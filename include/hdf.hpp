//
// Created by Jorrit Leenaarts on 2018-06-21.
//

#ifndef SHINY2_HDF_HPP
#define SHINY2_HDF_HPP

#include <string>
#include <vector>

#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"


using namespace std;

namespace hdf{
    /*
     * Class that handles reading (and soon writing) of hdf5 files
     */

    void init();
    void finalize();

    class Hdf{

    protected:

        hid_t fileID;
        hid_t faplID;
        int comm;
        bool isParallel = false;
        string filename;

    public:

        // serial open
        void open(string fname, string mode);
        // parallel open
        void open(string fname, string mode, MPI_Comm communicator);
        void close();

        bool isDatasetPresent(string name);
        // serial read (each process reads the same data)
        void readParallel(const string& name, const string& dtype, void* theData);
        // parallel read (each process reads the subdomain indicated by start and block
        void readParallel(const string& name, const string& dtype, void* theData,
                          const vector<int>& start, const vector<int>& block);

        // *********************************************************************

    };

}

#endif //SHINY2_HDF_HPP
