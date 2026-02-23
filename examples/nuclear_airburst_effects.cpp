/**
 * @file nuclear_airburst_effects.cpp
 * @brief Standalone driver for the NuclearAirburstEffects calculator
 *
 * Reads a config file, computes all airburst effects, and writes
 * summary / CSV / grid outputs.  No PETSc dependency — can run without MPI.
 *
 * Usage:
 *   ./nuclear_airburst_effects -c config/nuclear_airburst_berkeley_100kt.config
 *   ./nuclear_airburst_effects --generate-template my_scenario.config
 *
 * @author FSRM Development Team
 */

#include "NuclearAirburstEffects.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

static void printUsage(const char* prog)
{
  std::cout
    << "Usage:\n"
    << "  " << prog << " -c <config_file>             Run scenario from config\n"
    << "  " << prog << " --generate-template <file>    Write default template config\n"
    << "  " << prog << " --help                        Show this message\n\n"
    << "Example:\n"
    << "  " << prog << " -c config/nuclear_airburst_berkeley_100kt.config\n\n";
}

// Minimal INI reader for the [OUTPUT] section
struct OutputConfig
{
  bool   radial_profile      = true;
  std::string radial_file    = "radial_profile.csv";
  double radial_r_min        = 100.0;
  double radial_r_max        = 50000.0;
  int    radial_n            = 500;

  bool   grid_csv            = true;
  std::string grid_file      = "grid_output.csv";
  double grid_half_extent    = 25000.0;
  int    grid_n              = 200;

  bool   print_summary       = true;
};

static std::string trim(const std::string& s)
{
  auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

static OutputConfig parseOutputSection(const std::string& path)
{
  OutputConfig out;
  std::ifstream f(path);
  if (!f.is_open()) return out;

  std::string section, line;
  while (std::getline(f, line))
  {
    line = trim(line);
    if (line.empty() || line[0] == '#' || line[0] == ';') continue;
    if (line.front() == '[' && line.back() == ']')
    {
      section = trim(line.substr(1, line.size() - 2));
      continue;
    }
    if (section != "OUTPUT") continue;

    auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    std::string key = trim(line.substr(0, eq));
    std::string val = trim(line.substr(eq + 1));
    auto cp = val.find('#');
    if (cp != std::string::npos) val = trim(val.substr(0, cp));

    try
    {
      if      (key == "radial_profile")      out.radial_profile   = (val == "true" || val == "1");
      else if (key == "radial_profile_file") out.radial_file      = val;
      else if (key == "radial_r_min_m")      out.radial_r_min     = std::stod(val);
      else if (key == "radial_r_max_m")      out.radial_r_max     = std::stod(val);
      else if (key == "radial_n_points")     out.radial_n         = std::stoi(val);
      else if (key == "grid_csv")            out.grid_csv         = (val == "true" || val == "1");
      else if (key == "grid_csv_file")       out.grid_file        = val;
      else if (key == "grid_half_extent_m")  out.grid_half_extent = std::stod(val);
      else if (key == "grid_n_points")       out.grid_n           = std::stoi(val);
      else if (key == "print_summary")       out.print_summary    = (val == "true" || val == "1");
    }
    catch (...) {}
  }
  return out;
}

int main(int argc, char** argv)
{
  // --help
  if (argc < 2
      || std::strcmp(argv[1], "--help") == 0
      || std::strcmp(argv[1], "-h") == 0)
  {
    printUsage(argv[0]);
    return 0;
  }

  // --generate-template <file>
  if (std::strcmp(argv[1], "--generate-template") == 0)
  {
    std::string out_path = (argc >= 3) ? argv[2] : "airburst_template.config";
    FSRM::NuclearAirburstEffects calc;
    calc.writeTemplateConfig(out_path);
    return 0;
  }

  // -c <config>
  std::string config_path;
  for (int i = 1; i < argc; ++i)
  {
    if ((std::strcmp(argv[i], "-c") == 0 || std::strcmp(argv[i], "--config") == 0)
        && i + 1 < argc)
    {
      config_path = argv[i + 1];
      break;
    }
  }

  if (config_path.empty())
  {
    std::cerr << "Error: no config file specified.  Use -c <file>.\n";
    printUsage(argv[0]);
    return 1;
  }

  // ── Load scenario ─────────────────────────────────────────────────────────
  FSRM::NuclearAirburstEffects calc(config_path);
  OutputConfig out = parseOutputSection(config_path);

  std::cout << std::string(72, '=') << "\n"
            << "  FSRM — Nuclear Airburst Effects Calculator\n"
            << std::string(72, '=') << "\n\n"
            << "  Config:  " << config_path << "\n"
            << "  Yield:   " << calc.parameters().yield_kt << " kt\n"
            << "  HOB:     " << calc.parameters().burst_height_m << " m\n"
            << "  Location: " << calc.parameters().latitude << " N, "
                              << calc.parameters().longitude << " E  ("
                              << calc.parameters().location_name << ")\n\n";

  // ── Summary ───────────────────────────────────────────────────────────────
  if (out.print_summary)
  {
    calc.printSummary();
    std::cout << "\n";
  }

  // ── Radial profile ────────────────────────────────────────────────────────
  if (out.radial_profile)
  {
    calc.writeRadialProfile(out.radial_file, out.radial_r_min,
                            out.radial_r_max, out.radial_n);
  }

  // ── 2-D grid ──────────────────────────────────────────────────────────────
  if (out.grid_csv)
  {
    calc.writeGridCSV(out.grid_file, out.grid_half_extent, out.grid_n);
  }

  std::cout << "\nDone.\n";
  return 0;
}
