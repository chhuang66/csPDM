find_package(HDF5)
find_package(BLAS)
find_package(LAPACK)

if(HDF5_FOUND)
    set(HDF5_INCLUDE_DIR "/home/justin/anaconda2/include/")
	include_directories(${HDF5_INCLUDE_DIR})
	add_definitions(-DHAVE_HDF5)
endif()
set(libcstsrc sim/h5_file.cc prm/hdf5.cc prm/argu.cc prm/prm.cc arg.cc val.cc  sim/mt19937.cc sim/runner.cc)
add_library(cst ${libcstsrc})
target_link_libraries(cst ${HDF5_LIBRARIES})
set(libpdmsrc pdm/dgm.cpp pdm/Matrices.cpp pdm/nrn_model.cpp pdm/pdm.cpp pdm/filt.cpp)
add_library(pdm ${libpdmsrc})
target_link_libraries(pdm ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

