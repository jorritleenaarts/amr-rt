//
// Created by Jorrit Leenaarts on 2018-06-21.
//
#include <string>

#include <iostream>
#include "ftools.hpp"
#include "hdf.hpp"

using namespace std;

// *********************************************************************

void hdf::init(){
    if (H5open()){
        std::cout << "Error initializing HDF5 using H5open. Stopping execution." << std::endl;
        exit(1);
    }
}

// *********************************************************************

void hdf::finalize() {
    if (H5close()){
        std::cout << "Error closing HDF5 using H5close." << std::endl;
    }
}

// *********************************************************************

void hdf::Hdf::open(const string fname, const string mode) {
/*
 * open hdf5 file for single process access , mode is 'r','rw', or 'n'
 */

    filename = fname;

    if (!ftl::file_exists(fname) && (mode == "r" || mode == "rw")){
        ftl::message("hdf::Hdf::open: file " + fname
                     + " does not exists and access request is r or rw", 2);
    }

    const char* cfname = fname.c_str();
    if (mode == "r"){
        fileID = H5Fopen(cfname, H5F_ACC_RDONLY, H5P_DEFAULT);
    } else if (mode == "rw") {
        fileID = H5Fopen(cfname, H5F_ACC_RDWR, H5P_DEFAULT);
    } else if (mode == "n") {
        fileID = H5Fcreate (cfname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

}

// *********************************************************************

void hdf::Hdf::open(string fname, string mode, MPI_Comm communicator) {
/*
 * open hdf5 file for parallel access , mode is 'r','rw', or 'n'
 */

    if (!ftl::file_exists(fname) && (mode == "r" || mode == "rw")){
        ftl::message("hdf::Hdf::open: file " + fname
                     + " does not exists and access request is r or rw", 2);
    }

    hid_t faplID = H5Pcreate(H5P_FILE_ACCESS);
    MPI_Info info = MPI_INFO_NULL;
    H5Pset_fapl_mpio(faplID, communicator, info);

    const char* filename = fname.c_str();
    if (mode == "r"){
        fileID = H5Fopen(filename, H5F_ACC_RDONLY, faplID);
    } else if (mode == "rw") {
        fileID = H5Fopen(filename, H5F_ACC_RDWR, faplID);
    } else if (mode == "n") {
        fileID = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, faplID);
    } else {
        ftl::message("hdf::Hdf::open: mode '"+ mode + "' is not supported",2);
    }

    H5Pclose(faplID);
    isParallel = true;

}



void hdf::Hdf::close() {
    if (H5Fclose(fileID) < 0){
        ftl::message("hdf::Hdf::close: A problem occured when closing file  "
                     + filename, 2);
    }
}

// *********************************************************************

bool hdf::Hdf::isDatasetPresent(string name) {
/*
 * Determine whether a dataset exists at the top level  of the file
 */

    if (H5LTfind_dataset( fileID, name.c_str())){
        return true;
    } else {
        return false;
    }

}

// ********************************************************************

void hdf::Hdf::readParallel(const string& name, const string& dtype, void* theData,
                            const vector<int>& start2, const vector<int>& block2) {
/*
 * reads a 3D or 4D array, each process reads data for its own spatial/frequency subdomain
 */

    int s = start2.size();
    int s2 = block2.size();

    if (s != s2) {
        ftl::message("hdf::Hdf::readParallel: start2 and block2 are not the same size", 2);
    }

    hsize_t count[s];
    hsize_t stride[s];
    hsize_t start[s];
    hsize_t block[s];

    for (int i=0; i<s; i++){
        count[i] = 1;
        stride[i] = 1;
        start[i] = start2[i];
        block[i] = block2[i];
    }

    hid_t memtype_id;
    if (dtype == "int"){
        memtype_id = H5T_NATIVE_INT;
    } else if (dtype == "float"){
        memtype_id = H5T_NATIVE_FLOAT;
    } else if (dtype == "double"){
        memtype_id = H5T_NATIVE_DOUBLE;
    } else {
        ftl::message("hdf::Hdf::readParallel, datatype " + dtype + " not supported", 2);
    }

    hid_t dataset_id = H5Dopen(fileID, name.c_str(), H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, stride, count, block);
    hid_t memspace_id = H5Screate_simple(s, block, NULL);
    herr_t err = H5Dread( dataset_id, memtype_id, memspace_id, dataspace_id, H5P_DEFAULT, theData);
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

}

// ********************************************************************

void hdf::Hdf::readParallel(const string& name, const string& dtype, void* theData){
/*
 * read 1D array, all processes read all data in case of parallel IO
 */

    hid_t dataset_id = H5Dopen(fileID, name.c_str(), H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    hsize_t n = H5Sget_simple_extent_npoints(dataspace_id);

    hid_t memtype_id;
    if (dtype == "int"){
        memtype_id = H5T_NATIVE_INT;
    } else if (dtype == "float"){
        memtype_id = H5T_NATIVE_FLOAT;
    } else if (dtype == "double"){
        memtype_id = H5T_NATIVE_DOUBLE;
    } else {
        ftl::message("hdf::Hdf::readParallel, datatype " + dtype + " not supported", 2);
    }

    const int rank = 1;
    hid_t memspace_id = H5Screate_simple(rank, &n, NULL);
    herr_t err = H5Dread( dataset_id, memtype_id, memspace_id, dataspace_id, H5P_DEFAULT, theData);

    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

}

// ********************************************************************
