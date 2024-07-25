#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module is a collection of classes that provide a
friendlier interface to MPI (through mpi4py). They help
allocate local tasks/data from global tasks/data and gather
global data (from all processors).

Although general, this module was only implemented to
work with the convergence evaluation framework. More work
is needed to make this appropriate for general use.
"""
# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

from collections import OrderedDict
import importlib


class MPIInterface:
    __have_mpi__ = None

    def __init__(self):
        if MPIInterface.__have_mpi__ is None:
            # This is trying to import mpy4py.MPI, and setting a flag to indicate
            # if it succeeds or not.
            # we do this here instead of at the module level, because we only want
            # to do the import if an MPIInterface is ever requested.
            try:
                # try the import (the 'globals()' makes sure it is imported
                # in the module space and not local to the __init__ method)
                globals()["MPI"] = importlib.import_module("mpi4py.MPI")
                # import succeeded
                MPIInterface.__have_mpi__ = True
            except ImportError:
                # import failed (e.g., no mpi4py installed)
                MPIInterface.__have_mpi__ = False

        self._comm = None
        self._size = None
        self._rank = None

        if self.have_mpi:
            self._comm = MPI.COMM_WORLD  # pylint: disable=undefined-variable
            self._size = self._comm.Get_size()
            self._rank = self._comm.Get_rank()

    @property
    def have_mpi(self):
        assert MPIInterface.__have_mpi__ is not None
        return MPIInterface.__have_mpi__

    @property
    def comm(self):
        return self._comm

    @property
    def rank(self):
        return self._rank

    @property
    def size(self):
        return self._size


class ParallelTaskManager:
    def __init__(self, n_total_tasks, mpi_interface=None):
        if mpi_interface is None:
            self._mpi_interface = MPIInterface()
        else:
            self._mpi_interface = mpi_interface
        self._n_total_tasks = n_total_tasks

        if not self._mpi_interface.have_mpi:
            self._local_map = range(n_total_tasks)
        else:
            size = self._mpi_interface.size

            # there must be a better way to do this
            # find which entries in global correspond
            # to this process (want them to be contiguous
            # for the MPI Allgather calls later
            local_N = [0 for i in range(self._mpi_interface.size)]
            for i in range(n_total_tasks):
                process_i = i % size
                local_N[process_i] += 1

            start = 0
            end = None
            for i, v in enumerate(local_N):
                if i == self._mpi_interface.rank:
                    end = start + v
                    break
                else:
                    start += v

            self._local_map = list(range(start, end))

    def is_root(self):
        if not self._mpi_interface.have_mpi or self._mpi_interface.rank == 0:
            return True
        return False

    # ToDo: fix the parallel task manager to handle dictionaries as well as lists
    def global_to_local_data(self, global_data):
        if isinstance(global_data, list):
            local_data = list()
            assert len(global_data) == self._n_total_tasks
            for i in self._local_map:
                local_data.append(global_data[i])
            return local_data
        elif isinstance(global_data, OrderedDict):
            local_data = OrderedDict()
            assert len(global_data) == self._n_total_tasks
            idx = 0
            for k, v in global_data.items():
                if idx in self._local_map:
                    local_data[k] = v
                idx += idx
            return local_data
        raise ValueError(
            "Unknown type passed to global_to_local_data. Expected list or OrderedDict."
        )

    def allgather_global_data(self, local_data):
        assert len(local_data) == len(self._local_map)
        if not self._mpi_interface.have_mpi:
            return list(local_data)

        # PYLINT-TODO-FIX fix the error due to the
        # non-existing global_data_list_of_lists variable
        # pylint: disable=undefined-variable
        return self._stack_global_data(global_data_list_of_lists)

    def gather_global_data(self, local_data):
        assert len(local_data) == len(self._local_map)
        if not self._mpi_interface.have_mpi:
            return list(local_data)

        comm = self._mpi_interface.comm
        global_data_list_of_lists = comm.gather(local_data)

        if global_data_list_of_lists is not None:
            return self._stack_global_data(global_data_list_of_lists)

        assert self.is_root() is False
        return None

    def _stack_global_data(self, global_data_list_of_lists):
        # stack the list of lists into one global data list
        # ToDo: test that this is equivalent to [d for sublist in global_data_list_of_lists for d in sublist]
        global_data = list()
        for i in range(self._mpi_interface.size):
            global_data.extend(global_data_list_of_lists[i])
        return global_data
