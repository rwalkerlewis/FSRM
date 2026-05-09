/**
 * @file RadialLagrangianOutput.cpp
 * @brief HDF5 / XDMF spatial-profile writer implementation.
 *
 * Uses the HDF5 C API directly. The CMake configuration links HDF5
 * unconditionally when the library is found; this file falls back to
 * a no-op (returning true with a stderr warning) when HDF5 is absent
 * so the build does not regress in HDF5-less environments.
 */

#include "domain/explosion/RadialLagrangianOutput.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#if defined(__has_include)
#  if __has_include(<hdf5.h>)
#    include <hdf5.h>
#    define FSRM_HAVE_HDF5_H 1
#  endif
#endif

namespace FSRM {

#ifdef FSRM_HAVE_HDF5_H

namespace
{

bool writeDataset1D(hid_t group, const std::string& name,
                    const std::vector<double>& data)
{
    if (data.empty()) return true;
    hsize_t dims[1] = {static_cast<hsize_t>(data.size())};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) return false;
    hid_t dset = H5Dcreate2(group, name.c_str(), H5T_NATIVE_DOUBLE,
                            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) { H5Sclose(space); return false; }
    herr_t status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL,
                             H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(dset);
    H5Sclose(space);
    return status >= 0;
}

} // namespace

bool writeRadialProfilesHDF5(
    const std::string& path,
    const std::vector<RadialLagrangianSolver::RadialProfile>& profiles)
{
    if (profiles.empty()) return true;

    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                           H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "writeRadialProfilesHDF5: cannot create " << path
                  << "\n";
        return false;
    }

    // /time
    {
        std::vector<double> times;
        times.reserve(profiles.size());
        for (const auto& p : profiles) times.push_back(p.time);
        writeDataset1D(file, "time", times);
    }

    // /num_cells (scalar). Use the first profile's N.
    {
        const int N = static_cast<int>(profiles.front().rho.size());
        hid_t scalar = H5Screate(H5S_SCALAR);
        hid_t dset = H5Dcreate2(file, "num_cells", H5T_NATIVE_INT,
                                scalar, H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT);
        if (dset >= 0) {
            H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, &N);
            H5Dclose(dset);
        }
        H5Sclose(scalar);
    }

    // /profiles/<i>/...
    hid_t profiles_group = H5Gcreate2(file, "profiles", H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
    for (size_t i = 0; i < profiles.size(); ++i) {
        std::ostringstream name;
        name << std::setfill('0') << std::setw(6) << i;
        hid_t g = H5Gcreate2(profiles_group, name.str().c_str(),
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (g < 0) continue;
        const auto& p = profiles[i];
        writeDataset1D(g, "r", p.r_face);
        writeDataset1D(g, "r_cell", p.r_cell);
        writeDataset1D(g, "v_r", p.v_r);
        writeDataset1D(g, "rho", p.rho);
        writeDataset1D(g, "p", p.p);
        writeDataset1D(g, "sigma_rr", p.sigma_rr);
        writeDataset1D(g, "sigma_tt", p.sigma_tt);
        writeDataset1D(g, "eps_p", p.eps_p);
        writeDataset1D(g, "damage", p.damage);
        writeDataset1D(g, "yield_indicator", p.yield_indicator);
        H5Gclose(g);
    }
    H5Gclose(profiles_group);
    H5Fclose(file);
    return true;
}

#else // !FSRM_HAVE_HDF5_H

bool writeRadialProfilesHDF5(
    const std::string& path,
    const std::vector<RadialLagrangianSolver::RadialProfile>& profiles)
{
    if (!profiles.empty()) {
        std::cerr << "writeRadialProfilesHDF5: HDF5 not compiled in; "
                     "skipping " << path << "\n";
    }
    return true;
}

#endif // FSRM_HAVE_HDF5_H

bool writeRadialProfilesXDMF(
    const std::string& xdmf_path,
    const std::string& h5_basename,
    const std::vector<RadialLagrangianSolver::RadialProfile>& profiles)
{
    if (profiles.empty()) return true;

    std::ofstream xml(xdmf_path);
    if (!xml) {
        std::cerr << "writeRadialProfilesXDMF: cannot create "
                  << xdmf_path << "\n";
        return false;
    }

    const int N = static_cast<int>(profiles.front().rho.size());
    const int Nface = N + 1;

    xml << "<?xml version=\"1.0\" ?>\n"
        << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
        << "<Xdmf Version=\"3.0\">\n"
        << "  <Domain>\n"
        << "    <Grid Name=\"NearFieldProfileSeries\" "
           "GridType=\"Collection\" CollectionType=\"Temporal\">\n";

    for (size_t i = 0; i < profiles.size(); ++i) {
        std::ostringstream gname;
        gname << std::setfill('0') << std::setw(6) << i;
        const std::string g = gname.str();

        xml << "      <Grid Name=\"profile_" << g
            << "\" GridType=\"Uniform\">\n"
            << "        <Time Value=\"" << profiles[i].time << "\"/>\n"
            << "        <Topology TopologyType=\"Polyline\" "
               "NumberOfElements=\"" << N << "\" "
               "NodesPerElement=\"2\">\n"
            << "          <DataItem Dimensions=\"" << N << " 2\" "
               "Format=\"XML\" DataType=\"Int\">\n";
        for (int j = 0; j < N; ++j) {
            xml << "            " << j << " " << (j + 1) << "\n";
        }
        xml << "          </DataItem>\n"
            << "        </Topology>\n"
            << "        <Geometry GeometryType=\"X_Y_Z\">\n"
            << "          <DataItem Dimensions=\"" << Nface
            << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"
            << "            " << h5_basename << ":/profiles/" << g << "/r\n"
            << "          </DataItem>\n"
            << "          <DataItem Dimensions=\"" << Nface
            << "\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n";
        for (int j = 0; j < Nface; ++j) xml << "            0.0\n";
        xml << "          </DataItem>\n"
            << "          <DataItem Dimensions=\"" << Nface
            << "\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n";
        for (int j = 0; j < Nface; ++j) xml << "            0.0\n";
        xml << "          </DataItem>\n"
            << "        </Geometry>\n";

        // Cell-centred attributes: pressure, sigma_rr, sigma_tt, eps_p,
        // damage, yield_indicator, density.
        auto attr = [&](const std::string& name, const std::string& dset) {
            xml << "        <Attribute Name=\"" << name
                << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
                << "          <DataItem Dimensions=\"" << N
                << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"
                << "            " << h5_basename << ":/profiles/" << g
                << "/" << dset << "\n"
                << "          </DataItem>\n"
                << "        </Attribute>\n";
        };
        attr("rho", "rho");
        attr("p", "p");
        attr("sigma_rr", "sigma_rr");
        attr("sigma_tt", "sigma_tt");
        attr("eps_p", "eps_p");
        attr("damage", "damage");
        attr("yield_indicator", "yield_indicator");

        xml << "      </Grid>\n";
    }

    xml << "    </Grid>\n"
        << "  </Domain>\n"
        << "</Xdmf>\n";

    xml.close();
    return true;
}

} // namespace FSRM
