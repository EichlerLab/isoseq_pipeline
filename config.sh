module load mpc/0.8.2
module load mpfr/3.1.0
module load gmp/5.0.2
module load gcc/4.9.1
module load cmake/3.5.1
module load boost/1.58.0
module load zlib/1.2.6
module load libtool/2.4
module load autoconf/2.69
module load automake/1.14
module load curl/7.48.0
module load bamtools/2.3.0
module load samtools/0.1.19
module load perl/5.14.2

export BOOST_ROOT=$MOD_GSBOOST_DIR

export LD_LIBRARY_PATH=${TOP}/blasr_install/blasr/libcpp/alignment:${TOP}/blasr_install/blasr/libcpp/hdf:${TOP}/blasr_install/blasr/libcpp/pbdata:${TOP}/blasr_install/hdf5/hdf5-1.8.16-linux-centos6-x86_64-gcc447-shared/lib/:$LD_LIBRARY_PATH

bash --init-file /net/eichler/vol18/zevk/great_apes/iso_seq/cc2_analysis/pitchfork/deployment/setup-env.sh
