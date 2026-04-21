/* ============================================================
   jpl_fetch.cpp  --  SolarHPC
   Reads data/initial_conditions.dat into BodyState structs
   and provides flat C-arrays for Fortran calls.
   ============================================================ */
#include "solar_types.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>

/* ── load_initial_conditions ─────────────────────────────── */
bool load_initial_conditions(const std::string& filepath,
                              std::vector<BodyState>& bodies)
{
    std::ifstream f(filepath);
    if (!f.is_open()) {
        std::cerr << "  [ERROR] Cannot open: " << filepath << "\n";
        std::cerr << "  Run: python3 fetch_jpl.py\n";
        return false;
    }

    bodies.clear();
    std::string line;
    int line_no = 0;

    while (std::getline(f, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        BodyState b{};
        std::string name;
        if (!(ss >> name >> b.mass_kg
                 >> b.pos[0] >> b.pos[1] >> b.pos[2]
                 >> b.vel[0] >> b.vel[1] >> b.vel[2])) {
            std::cerr << "  [ERROR] Parse failed at line " << line_no << "\n";
            return false;
        }
        strncpy(b.name, name.c_str(), 31);
        b.name[31]    = '\0';
        b.body_index  = (int)bodies.size();
        bodies.push_back(b);
    }

    if ((int)bodies.size() < N_BODIES_CPP) {
        std::cerr << "  [ERROR] Only " << bodies.size()
                  << " bodies parsed (expected " << N_BODIES_CPP << ")\n";
        return false;
    }
    return true;
}

/* ── flatten_to_fortran ──────────────────────────────────── */
/* Fortran arrays are column-major: pos(3,N)
   pos_out[0 + i*3] = bodies[i].pos[0]  (x of body i)
   pos_out[1 + i*3] = bodies[i].pos[1]  (y)
   pos_out[2 + i*3] = bodies[i].pos[2]  (z)
   Fortran body index i (1-based) maps to C++ i-1.            */
void flatten_to_fortran(const std::vector<BodyState>& bodies,
                         double* pos_out,
                         double* vel_out,
                         double* mass_out)
{
    for (int i = 0; i < (int)bodies.size(); ++i) {
        pos_out[0 + i*3] = bodies[i].pos[0];
        pos_out[1 + i*3] = bodies[i].pos[1];
        pos_out[2 + i*3] = bodies[i].pos[2];
        vel_out[0 + i*3] = bodies[i].vel[0];
        vel_out[1 + i*3] = bodies[i].vel[1];
        vel_out[2 + i*3] = bodies[i].vel[2];
        mass_out[i]       = bodies[i].mass_kg;
    }
}

/* ── print_body_table ────────────────────────────────────── */
void print_body_table(const std::vector<BodyState>& bodies)
{
    const char* sep = "  +-----------+--------------+-------------+-------------+-------------+";
    std::cout << sep << "\n";
    std::cout << "  | Body      | Mass (kg)    | x (AU)      | y (AU)      | z (AU)      |\n";
    std::cout << sep << "\n";
    for (const auto& b : bodies) {
        std::cout << std::fixed << std::setprecision(4)
                  << "  | " << std::left  << std::setw(10) << b.name
                  << "| " << std::right << std::setw(12) << std::scientific
                            << std::setprecision(3) << b.mass_kg
                  << " | " << std::setw(11) << std::fixed << std::setprecision(6)
                            << b.pos[0]
                  << " | " << std::setw(11) << b.pos[1]
                  << " | " << std::setw(11) << b.pos[2]
                  << " |\n";
    }
    std::cout << sep << "\n\n";
}

/* ── bodies_initial_energy ───────────────────────────────── */
double bodies_initial_energy(const std::vector<BodyState>& bodies)
{
    double KE = 0.0, PE = 0.0;
    for (const auto& b : bodies) {
        double v2 = 0.0;
        for (int k = 0; k < 3; ++k)
            v2 += b.vel[k] * b.vel[k];
        KE += 0.5 * b.mass_kg * v2 * AU_PER_DAY_TO_M_PER_S * AU_PER_DAY_TO_M_PER_S;
    }
    for (int i = 0; i < (int)bodies.size()-1; ++i) {
        for (int j = i+1; j < (int)bodies.size(); ++j) {
            double dx = (bodies[j].pos[0]-bodies[i].pos[0]) * AU_TO_M;
            double dy = (bodies[j].pos[1]-bodies[i].pos[1]) * AU_TO_M;
            double dz = (bodies[j].pos[2]-bodies[i].pos[2]) * AU_TO_M;
            double r  = std::sqrt(dx*dx+dy*dy+dz*dz);
            PE -= G_SI * bodies[i].mass_kg * bodies[j].mass_kg / r;
        }
    }
    return KE + PE;
}
