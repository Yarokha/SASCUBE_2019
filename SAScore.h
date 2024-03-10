#pragma once
#pragma once
#ifndef SAScore_H_
#define SAScore_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <fstream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <ctime>
#include <valarray>
#include <map>
#include <boost/math/special_functions/sinc.hpp>
#include <PDBParser.h>

#define knum_threads std::thread::hardware_concurrency()
#define boost_sinc_pi boost::math::sinc_pi 

using namespace std::complex_literals;

class SAScore
{
private:
    //threads
//    int knum_threads = std::thread::hardware_concurrency()  ;
    std::thread *thread_ = new std::thread[knum_threads];
    std::mutex mtx_;


    const std::map <double, double> kmolcoef_= {
        {0.02, 0.052},
        {0.03, 0.120}
   };
    const std::map <double, double> khydcoef_ = {
       {0.02, 0.05},
       {0.03, 0.115}
    };
    
    std::string path_folder_="";
    //structure
    //amount of components for initialization
    enum class AmountCubeComponents_ {
        kCoord = 3,
        kCoordState = 4,
    };

    //Number of component of cube
    enum class CubeComponents_ {
        kX = 0,
        kY = 1,
        kZ = 2,
        kState = 3,
    };

    //Number of component of pdb file
    enum class PDBComponents_ {
        kX = 0,
        kY = 1,
        kZ = 2,
        kRadius = 3,
        kCharge = 4,
        ka0 = 5,
        ka5 = 10,
    };

    //states of cube
    enum struct CubeState_ {
        kEmpty = 0,
        kFilled = 1,
        kBoundary = 2,
    };

    //template of atom
    struct AtomTemplate_ {
        double r;
        std::vector< std::vector<int> > grid_atom = std::vector< std::vector<int> >(0);
    };
    struct Parallelepiped_ {
        std::vector<int> first_cube = std::vector<int>(0);
        std::vector<int> last_cube = std::vector<int>(0);
        std::vector<double> coordinates = std::vector<double>(int(AmountCubeComponents_::kCoord));
        int length = 0;
        double x_length = 0, y_length = 0, z_length = 0, volume = 0;
    }big_piped_;

//    struct PARALLELEPIPED {
//            std::vector<double>
//                x = std::vector<double>(2),
//                y = std::vector<double>(2),
//                z = std::vector<double>(2);
//    };

    struct QSpace_ {
        std::vector <std::vector <double>> q0_x = std::vector <std::vector <double>>(0);
        std::vector <std::vector <double>> q0_y = std::vector <std::vector <double>>(0);
        std::vector <std::vector <double>> q0_z = std::vector <std::vector <double>>(0);
        int N_z = 70, N_phi = 70;
    }q_space_;

    //variables
    const double step_;
    double chi_ = 0;
    std::complex <double> rho_0_ = 0;
    std::vector< std::vector<double> > pdb_data_;
    std::vector <AtomTemplate_> atom_templates_ = std::vector <AtomTemplate_>(0);
    struct ThreadMemory_ {
        std::vector< std::vector<double> > pdb_data = std::vector< std::vector<double> >(0);
        std::vector< std::vector<int> > hydrated_volume = std::vector< std::vector<int> >(0);
        std::vector< std::vector<int> > molecular_volume = std::vector< std::vector<int> >(0);

    };

    std::vector <ThreadMemory_> threads_memory_ = std::vector <ThreadMemory_>(knum_threads);
    std::vector <Parallelepiped_> parallelepipeds_ = std::vector <Parallelepiped_>(0);
    std::vector< std::vector<int> > total_hydrated_volume_ = std::vector< std::vector<int> >(0);
    std::vector< std::vector<int> >total_molecular_volume_ = std::vector< std::vector<int> >(0);
    std::vector< std::vector <double > > I_ = std::vector< std::vector <double> >(0);

    //functions
    //for sorting
    static bool XYZ_(const std::vector <int>& vec1, const std::vector <int>& vec2);
    static bool YZX_(const std::vector <int>& vec1, const std::vector <int>& vec2);
    static bool ZXY_(const std::vector <int>& vec1, const std::vector <int>& vec2);
    static bool XYZState_(const std::vector <int>& vec1, const std::vector <int>& vec2);
    static bool StateInv_(const std::vector <int>& vec1, const std::vector <int>& vec2);
    static bool dXYZ_(const std::vector <double>& vec1, const std::vector <double>& vec2);
    static bool ByQ_(const std::vector <double>& vec1, const std::vector <double>& vec2);
    //for determining state of cube
    bool DetermineCubeStateXYZ_(const std::vector<int>& cube1, const std::vector<int>& cube2);
    bool DetermineCubeStateYZX_(const std::vector<int>& cube1, const std::vector<int>& cube2);
    bool DetermineCubeStateZXY_(const std::vector<int>& cube1, const std::vector<int>& cube2);
	
	template<typename T>
    void WriteInFile_(std::vector< std::vector<T> >& data, std::string filename);

    bool RadiusCheck_(const double& r);
    void AtomGridding_(const double& r);
    void AtomGridding_(const double& r, const double& delta);
    void AtomTemplatesBuilder_(std::vector< std::vector<double> >& data);
    void NewSpaceBuilder_(const std::vector<double>& data, std::vector< std::vector<int> >& volume);
    AtomTemplate_* AtomTemplateSearcher_(const double& r);
    void HydratedVolumeBuilder_(std::vector< std::vector<double> >& data);
    bool CoordCompare_(const std::vector<int>& a, const std::vector<int>& b);
    void DeleteEqualCoord_(std::vector< std::vector<int> >& data);
    void BoundaryBuilder_(std::vector< std::vector<int> >& data);
    void MolecularVolumeBuilder_(std::vector< std::vector<int> >& data);
    void SolventTemplateApply_(const std::vector<int>& data, std::vector<std::vector<int> >& output);
    void SaveUniqCoord_(std::vector< std::vector<int> >& data);
    void ParallelepipedSeacher_(std::vector< std::vector<int> >& data);
    void CalculatingQ0_();
    void IntensityCurveCalc_(double Min_q = 0.05, double Max_q = 9, int N_q = 100);
    void IntensityCurveCalc_(const std::vector<double>& q);
    void CalculatingChi_(std::vector< std::vector<double> >& data);
    double CalcFactor_(double& q, std::vector<double>& atom);
    int NumBorderCubes_(std::vector< std::vector<int> >& data);
    void BigPipedParam();
    void IntensityCurveCalcBigPiped_(double Min_q, double Max_q, int N_q);

public:
    SAScore(std::vector< std::vector<double> >& pdb_data, const double& solvent_radius, const double& step, const double& delta);
    ~SAScore();
    void CoreStart();
    void PipedLoad(std::string file_name);
    void CalcIntensity(double Min_q = 0.05, double Max_q = 10, int N_q = 100, std::complex <double> rho0 = 334.0);
    void CalcIntensityBigPiped(double Min_q, double Max_q, int N_q, std::complex <double> rho0);
    double GetChi();
    inline void ProjectFolder(std::string path){path_folder_ = path; /*freopen((path_folder_ + "/log.txt").c_str(), "w", stdout);*/}
	double ProteinVolume();
    double RadiusGuinier();
    void CalculatingIntChi(std::vector< std::vector<double> >& data, std::complex <double> rho0 = 334.0);
    void SetQSpace(int N_z, int N_phi);
    std::vector < std::vector<int> > ProteinModel();
    struct dimension{
        double hydrated = 0;
        double molecular = 0;
    } volume, square;

    PDBParser::PARALLELEPIPED big_piped;
};

#endif
