 #include "SAScore.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//sorting by coordinates in XYZ
bool SAScore::XYZ_(const std::vector <int>& vec1, const std::vector <int>& vec2) {
    for (int i = 0; i <= int(CubeComponents_::kZ); ++i) {
        if (vec1[i] == vec2[i]) continue;
        return vec1[i] < vec2[i];
    }
    return false;
}

bool SAScore::ByQ_(const std::vector <double>& vec1, const std::vector <double>& vec2) {
    return vec1[0] < vec2[0];
}

//sorting by coordinates in XYZ
bool SAScore::dXYZ_(const std::vector <double>& vec1, const std::vector <double>& vec2) {
    for (int i = 0; i <= int(CubeComponents_::kZ); ++i) {
        if (vec1[i] == vec2[i]) continue;
        return vec1[i] < vec2[i];
    }
    return false;
}

//sorting by coordinates in YZX
bool SAScore::YZX_(const std::vector <int>& vec1, const std::vector <int>& vec2) {
    for (int i = 1; i <= int(CubeComponents_::kZ); ++i) {
        if (vec1[i] == vec2[i]) continue;
        return vec1[i] < vec2[i];
    }
    return vec1[int(CubeComponents_::kX)] < vec2[int(CubeComponents_::kX)];
}

//sorting by coordinates in ZXY
bool SAScore::ZXY_(const std::vector <int>& vec1, const std::vector <int>& vec2) {
    if (vec1[int(CubeComponents_::kZ)] == vec2[int(CubeComponents_::kZ)]) {
        for (int i = 0; i <int(CubeComponents_::kZ); ++i) {
            if (vec1[i] == vec2[i]) continue;
            return vec1[i] < vec2[i];
        }
    }
    return vec1[int(CubeComponents_::kZ)] < vec2[int(CubeComponents_::kZ)];
}

//sorting by coordinates in XYZ and state
bool SAScore::XYZState_(const std::vector <int>& vec1, const std::vector <int>& vec2) {
    for (int i = 0; i <= int(CubeComponents_::kState); ++i) {
        if (vec1[i] == vec2[i]) continue;
        return vec1[i] < vec2[i];
    }
    return vec1[int(CubeComponents_::kState)] < vec2[int(CubeComponents_::kState)];
}

//inversed sorting only by state
bool SAScore::StateInv_(const std::vector <int>& vec1, const std::vector <int>& vec2) {
    return vec1[int(CubeComponents_::kState)] > vec2[int(CubeComponents_::kState)];
}
///////////////////////////////////////////////////////////////////////////////////////////////////

//Gridding of a atom
void SAScore::AtomGridding_(const double& r) {
    int radius_int = round(r / step_);
    int coord[int(AmountCubeComponents_::kCoord)] = { 0 };
    double r_sq = pow(r, 2), test_radius_sq = 0;
    AtomTemplate_ atom_template;
    atom_template.r = r;
    std::vector <int> cube(int(AmountCubeComponents_::kCoordState), 0);
    for (coord[int(CubeComponents_::kX)] = 0; coord[int(CubeComponents_::kX)] < 2 * radius_int + 1; ++coord[int(CubeComponents_::kX)]) {
        for (coord[int(CubeComponents_::kY)] = 0; coord[int(CubeComponents_::kY)] < 2 * radius_int + 1; ++coord[int(CubeComponents_::kY)]) {
            for (coord[int(CubeComponents_::kZ)] = 0; coord[int(CubeComponents_::kZ)] < 2 * radius_int + 1; ++coord[int(CubeComponents_::kZ)]) {
                //creating new coordinates and calculating cube distance from the centre of sphere
                for (int i = 0; i < int(AmountCubeComponents_::kCoord); ++i) {
                    cube[i] = coord[i] - radius_int;
                    test_radius_sq += pow(cube[i] * step_, 2);
                }
                if (test_radius_sq <= r_sq) {
                    cube[int(CubeComponents_::kState)] = int(CubeState_::kFilled);
                    atom_template.grid_atom.push_back(cube);
                }
                test_radius_sq = 0;
            }
        }
    }
    atom_templates_.push_back(atom_template);
}

void SAScore::AtomGridding_(const double& r, const double& delta) {
    int radius_int = round(r / step_);
    int coord[int(AmountCubeComponents_::kCoord)] = { 0 };
    double r_sq = pow(r-delta, 2), test_radius_sq = 0;
    AtomTemplate_ atom_template;
    atom_template.r = r;
    std::vector <int> cube(int(AmountCubeComponents_::kCoordState), 0);
    for (coord[int(CubeComponents_::kX)] = 0; coord[int(CubeComponents_::kX)] < 2 * radius_int + 1; ++coord[int(CubeComponents_::kX)]) {
        for (coord[int(CubeComponents_::kY)] = 0; coord[int(CubeComponents_::kY)] < 2 * radius_int + 1; ++coord[int(CubeComponents_::kY)]) {
            for (coord[int(CubeComponents_::kZ)] = 0; coord[int(CubeComponents_::kZ)] < 2 * radius_int + 1; ++coord[int(CubeComponents_::kZ)]) {
                //creating new coordinates and calculating cube distance from the centre of sphere
                for (int i = 0; i < int(AmountCubeComponents_::kCoord); ++i) {
                    cube[i] = coord[i] - radius_int;
                    test_radius_sq += pow(cube[i] * step_, 2);
                }
                if (test_radius_sq <= r_sq) {
                    cube[int(CubeComponents_::kState)] = int(CubeState_::kFilled);
                    atom_template.grid_atom.push_back(cube);
                }
                test_radius_sq = 0;
            }
        }
    }
    atom_templates_.push_back(atom_template);


}

//checking existing atom template
bool SAScore::RadiusCheck_(const double& r) {
    for (auto& iter : atom_templates_)
        if (r == iter.r) {
            return true;
        }
    return false;
}

void SAScore::AtomTemplatesBuilder_(std::vector< std::vector<double> >& data) {
    for (unsigned int i = 0; i < data.size(); ++i) {
        if (!RadiusCheck_(data[i][int(PDBComponents_::kRadius)]))
            AtomGridding_(data[i][int(PDBComponents_::kRadius)]);
    }
}

SAScore::AtomTemplate_* SAScore::AtomTemplateSearcher_(const double& r) {
    for (auto& iter : atom_templates_)
        if (r == iter.r) {
            return &iter;
        }
    return NULL;
}

//building new space
void SAScore::NewSpaceBuilder_(const std::vector<double>& data, std::vector< std::vector<int> >& volume) {
    AtomTemplate_ *existing_atom_template = NULL;
    existing_atom_template = AtomTemplateSearcher_(data[int(PDBComponents_::kRadius)]);
    std::vector <int> centre_coord(int(AmountCubeComponents_::kCoord), 0);
    std::vector <int> new_coord(int(AmountCubeComponents_::kCoordState), 0);
    for (int i = 0; i < int(AmountCubeComponents_::kCoord); ++i) {
        centre_coord[i] = round(data[i] / step_);
    }
    for (unsigned int i = 0; i < existing_atom_template->grid_atom.size(); ++i) {
        for (int j = 0; j < int(AmountCubeComponents_::kCoord); ++j) {
            new_coord[j] = centre_coord[j] + existing_atom_template->grid_atom[i][j];
        }
        new_coord[int(CubeComponents_::kState)] =
            existing_atom_template->grid_atom[i][int(CubeComponents_::kState)];
        volume.push_back(new_coord);
    }
}

//comparing of coordinates
bool SAScore::CoordCompare_(const std::vector<int>& a, const std::vector<int>& b) {
    for (int i = 0; i < int(AmountCubeComponents_::kCoord); ++i) {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

//Deleting equal coordinate(s)
void SAScore::DeleteEqualCoord_(std::vector< std::vector<int> >& data) {
    std::vector< std::vector<int> > temp(0);
    if (!CoordCompare_(data[0], data[1])) temp.push_back(data[0]);
    for (unsigned int i = 1; i < data.size(); ++i) {
        if (!CoordCompare_(data[i - 1], data[i]))
            temp.push_back(data[i]);
    }
    data = temp;
}

bool SAScore::DetermineCubeStateZXY_(const std::vector<int>& cube1, const std::vector<int>& cube2) {
    int d = 0, c = 0;
    for (int i = 0; i < int(AmountCubeComponents_::kCoord); ++i) {
        c = abs(cube1[i] - cube2[i]);
        if (i % 2 == 0 && c > 0)
            return true;
        d += c;
    }
    if (d > 1) return true;
    else return false;
}

bool SAScore::DetermineCubeStateYZX_(const std::vector<int>& cube1, const std::vector<int>& cube2) {
    int d = 0, c = 0;
    for (int i = 0; i <= int(CubeComponents_::kZ); ++i) {
        c = abs(cube1[i] - cube2[i]);
        if (i > int(CubeComponents_::kX) && c > 0)
            return true;
        d += c;
    }
    if (d > 1) return true;
    else return false;
}

bool SAScore::DetermineCubeStateXYZ_(const std::vector<int>& cube1, const std::vector<int>& cube2) {
    int d = 0;
    for (int i = 0; i < int(AmountCubeComponents_::kCoord); ++i) {
        d += abs(cube1[i] - cube2[i]);
        if (i <= int(CubeComponents_::kY) && d>0)
            return true;
    }
    if (d > 1) return true;
    else return false;
}

void SAScore::BoundaryBuilder_(std::vector< std::vector<int> >& data) {
    std::cout << "    In Z-direction";
    sort(data.begin(), data.end(), XYZ_);
    for (unsigned int i = 0; i < data.size() - 1; ++i) {
        if (DetermineCubeStateXYZ_(data[i], data[i + 1]))
            data[i][int(CubeComponents_::kState)] =
            data[i + 1][int(CubeComponents_::kState)] =
            int(CubeState_::kBoundary);
    }
    std::cout << "  Done!!!" << std::endl;

    std::cout << "    In Y-direction";
    sort(data.begin(), data.end(), ZXY_);
    for (unsigned int i = 0; i < data.size() - 1; ++i) {
        if (DetermineCubeStateZXY_(data[i], data[i + 1]))
            data[i][int(CubeComponents_::kState)] =
            data[i + 1][int(CubeComponents_::kState)] =
            int(CubeState_::kBoundary);
    }
    std::cout << "  Done!!!" << std::endl;

    std::cout << "    In X-direction";
    sort(data.begin(), data.end(), YZX_);
    for (unsigned int i = 0; i < data.size() - 1; ++i) {
        if (DetermineCubeStateYZX_(data[i], data[i + 1]))
            data[i][int(CubeComponents_::kState)] =
            data[i + 1][int(CubeComponents_::kState)] =
            int(CubeState_::kBoundary);
    }
    std::cout << "  Done!!!" << std::endl;
}

void SAScore::HydratedVolumeBuilder_(std::vector< std::vector<double> >& data) {
    std::time_t tStart = time(0);
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Building of Hydrated Volume" << std::endl;
    std::cout << "Building of Atom Templates";
    AtomTemplatesBuilder_(data);
    std::cout << "   Done!!!" << std::endl;
    std::cout << "Building of New Space" << std::endl;;
    if (knum_threads > 1)
        std::cout << "Prepearing for multithreading";
    sort(data.begin(), data.end(), dXYZ_);
    bool *thread_state = new bool[knum_threads] {false};
    int l = std::round(data.size() / knum_threads) + 1;
    for (unsigned int i = 0; i < knum_threads; ++i) {
        unsigned int begin = i*l, end = (i + 1)*l;
        if (end >= data.size()) {
            end = data.size();
        }
        threads_memory_[i].pdb_data.insert(
            std::end(threads_memory_[i].pdb_data),
            data.begin() + begin, data.begin() + end);
        thread_state[i] = true;
        if (end == data.size())
            break;
    }

    if (knum_threads > 1)
        std::cout << "   Done!!!" << std::endl;

    for (unsigned int j = 0; j < knum_threads; ++j) {
        if (thread_state[j] == true) {
            mtx_.lock();
            std::cout << "Thread number " << j << " starts" << std::endl;
            mtx_.unlock();
            thread_[j] = std::thread([this](std::vector< std::vector<double> >& data,
                const int num_thread) {
                for (auto iter : data) {
                    NewSpaceBuilder_(iter, threads_memory_[num_thread].hydrated_volume);
                }
                std::sort(threads_memory_[num_thread].hydrated_volume.begin(),
                    threads_memory_[num_thread].hydrated_volume.end(), XYZState_);
                DeleteEqualCoord_(threads_memory_[num_thread].hydrated_volume);
                mtx_.lock();
                std::cout << "Thread number " << num_thread << " finished" << std::endl;
                mtx_.unlock();
            }, std::ref(threads_memory_[j].pdb_data), j);
        }
    }

    for (unsigned int j = 0; j < knum_threads; ++j) {
        if (thread_state[j] == true)
            thread_[j].join();
    }

    for (unsigned int j = 0; j < knum_threads; ++j) {
        total_hydrated_volume_.insert(
            std::end(total_hydrated_volume_),
            std::begin(threads_memory_[j].hydrated_volume),
            std::end(threads_memory_[j].hydrated_volume));
        threads_memory_[j].hydrated_volume.clear();
    }

    std::cout << "Sorting of Coordinates";
    std::sort(total_hydrated_volume_.begin(), total_hydrated_volume_.end(), XYZState_);
    std::cout << "   Done!!!" << std::endl;
    std::cout << "Deleting of Equal Coordinates";
    DeleteEqualCoord_(total_hydrated_volume_);
    std::cout << "   Done!!!" << std::endl;
    std::cout << "Building of Hydrated Surface:" << std::endl;
    BoundaryBuilder_(total_hydrated_volume_);
    std::cout << "Writting in file";
    volume.hydrated = total_hydrated_volume_.size()* pow(step_,3);
    square.hydrated = NumBorderCubes_(total_hydrated_volume_);
    WriteInFile_(total_hydrated_volume_, "Hydrated");
    std::cout << "   Done!!!" << std::endl;
    std::cout << "Finished for " << std::setprecision(5)
        << (time(0) - tStart)<< " seconds" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}
template<typename T>
void SAScore::WriteInFile_(std::vector< std::vector<T> > &data, std::string filename) {
    std::ofstream file(path_folder_ +"/"+ filename + ".txt");
    for (auto & iter : data) {
        for (auto i : iter) {
            file << i << ' ';
        }
        file << std::endl;
    }
    file.close();
}

void SAScore::SolventTemplateApply_(const std::vector<int>& data, std::vector<std::vector<int> >& output) {
    std::vector <int> cube(int(AmountCubeComponents_::kCoordState), 0);
    for (auto i : atom_templates_[0].grid_atom) {
        for (int j = 0; j < int(AmountCubeComponents_::kCoord); ++j) {
            cube[j] = data[j] + i[j];
        }
        cube[int(CubeComponents_::kState)] = int(CubeState_::kEmpty);
        output.push_back(cube);
    }
}

void SAScore::SaveUniqCoord_(std::vector< std::vector<int> >& data) {
    sort(data.begin(), data.end(), XYZ_);
    size_t l = data.size() - 1;
    std::vector < std::vector<int> > temp(0);
    if (!CoordCompare_(data[0], data[1])) {
        if (data[0][int(CubeComponents_::kState)] == 1) {
            temp.push_back(data[0]);
        }
    }
    for (unsigned int i = 1; i < l - 1; i += 1) {
        if (!CoordCompare_(data[i], data[i + 1]) && !CoordCompare_(data[i], data[i - 1])) {
            if (data[i][int(CubeComponents_::kState)] == 1) {
                temp.push_back(data[i]);
            }
        }
    }
    if (!CoordCompare_(data[l - 1], data[l])) {
        if (data[l][int(CubeComponents_::kState)] == 1) {
            temp.push_back(data[l]);
        }
    }
    data = temp;
}

void SAScore::MolecularVolumeBuilder_(std::vector< std::vector<int> >& data) {
    std::time_t tStart = time(0);
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Building of Molecular Volume" << std::endl;
    if (knum_threads > 1)
        std::cout << "Prepearing for multithreading";
    sort(data.begin(), data.end(), XYZ_);
    int l = std::round(data.size() / knum_threads) + 1;
    bool *thread_state = new bool[knum_threads] {false};
    for (unsigned int i = 0; i < knum_threads; ++i) {
        unsigned int begin = i*l, end = (i + 1)*l;
        if (end >= data.size()) {
            end = data.size();
        }
        threads_memory_[i].hydrated_volume.insert(
            std::end(threads_memory_[i].hydrated_volume),
            data.begin() + begin, data.begin() + end);
        threads_memory_[i].molecular_volume = threads_memory_[i].hydrated_volume;
        thread_state[i] = true;
        if (end == data.size())
            break;
    }
    data.clear();
    if (knum_threads > 1)
        std::cout << "   Done!!!" << std::endl;

    for (unsigned int j = 0; j < knum_threads; ++j) {
        mtx_.lock();
        std::cout << "Thread number " << j << " starts" << std::endl;
        mtx_.unlock();
        thread_[j] = std::thread([this](std::vector< std::vector<int> >& data,
            const int num_thread) {
            sort(data.begin(), data.end(), StateInv_);
            for (auto& i : data) {
                if (i[int(CubeComponents_::kState)] != int(CubeState_::kBoundary))
                    break;
                SolventTemplateApply_(i, threads_memory_[num_thread].molecular_volume);
            }
            std::sort(threads_memory_[num_thread].molecular_volume.begin(),
                threads_memory_[num_thread].molecular_volume.end(), XYZState_);
            mtx_.lock();
            std::cout << "Thread number " << num_thread << " finished" << std::endl;
            mtx_.unlock();
        }, std::ref(threads_memory_[j].hydrated_volume), j);
    }
    for (unsigned int j = 0; j < knum_threads; ++j) {
        if (thread_state[j] == true)
            thread_[j].join();
    }

    for (unsigned int j = 0; j < knum_threads; ++j) {
        total_molecular_volume_.insert(
            std::end(total_molecular_volume_),
            std::begin(threads_memory_[j].molecular_volume),
            std::end(threads_memory_[j].molecular_volume));
    }
    threads_memory_.clear();

    std::cout << "Saving Uniq Coordinates";
    SaveUniqCoord_(total_molecular_volume_);
    std::cout << "   Done!!!" << std::endl;

    std::cout << "Building of Molecular Surface:" << std::endl;
    BoundaryBuilder_(total_molecular_volume_);

    volume.molecular = total_molecular_volume_.size()*pow(step_,3);
    square.molecular = NumBorderCubes_(total_molecular_volume_);
    std::cout << "Writting in file";
    WriteInFile_(total_molecular_volume_, "Molecular");
    std::cout << "   Done!!!" << std::endl;

    std::cout << "Finished for " << std::setprecision(5)
        << (time(0) - tStart) << " seconds" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}

void SAScore::ParallelepipedSeacher_(std::vector< std::vector<int> >& data) {
    std::time_t tStart = time(0);
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Searching of Paralellepipeds in the Molecular Volume" << std::endl;
    sort(data.begin(), data.end(), XYZ_);
    Parallelepiped_ temp;
    temp.first_cube = data[0];
    int l = 0;
    for (auto iter = data.begin() + 1; iter != data.end(); ++iter) {
        if (DetermineCubeStateXYZ_(*(iter - 1), *iter)) {
            temp.last_cube = *(iter - 1);
            temp.length = temp.last_cube[int(CubeComponents_::kZ)] -
                temp.first_cube[int(CubeComponents_::kZ)] + 1;
            l += temp.length;
            for (int i = 0; i <= int(CubeComponents_::kZ); ++i) {
                temp.coordinates[i] = (temp.last_cube[i] + temp.first_cube[i])*step_ / 2.0;
            }
            parallelepipeds_.push_back(temp);
            temp.first_cube = *iter;
        }
    }
    temp.last_cube = data[data.size() - 1];
    temp.length = temp.last_cube[int(CubeComponents_::kZ)] -
        temp.first_cube[int(CubeComponents_::kZ)] + 1;
    for (int i = 0; i <= int(CubeComponents_::kZ); ++i) {
        temp.coordinates[i] = (temp.last_cube[i] + temp.first_cube[i])*step_ / 2.0;
    }
    l += temp.length;
    parallelepipeds_.push_back(temp);

    data.clear();

    std::cout << "Writting in file";

    std::ofstream pararel;
    pararel.open(path_folder_ + "/piped.txt");

    for (auto i : parallelepipeds_) {
    for (int j = 0; j < int(AmountCubeComponents_::kCoord); ++j) {
    pararel << i.coordinates[j] << ' ';
    }
    pararel << i.length << std::endl;
    }
    pararel.close();

    std::cout << "   Done!!!" << std::endl;

    std::cout << "Finished for " << std::setprecision(5)
        <<(time(0) - tStart) << " seconds" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}

SAScore::SAScore(std::vector< std::vector<double> >& pdb_data, const double& solvent_radius, const double& step, const double& delta) :step_(step), pdb_data_(pdb_data)
{
    AtomGridding_(solvent_radius, delta);
    for(auto& i: pdb_data_){
        i[3]+=solvent_radius;
    }
}

SAScore::~SAScore()
{
}

void SAScore::CoreStart() {
    std::cout << "Number of threads: " << knum_threads << std::endl;
    std::time_t tStart = time(0);
    HydratedVolumeBuilder_(pdb_data_);
    MolecularVolumeBuilder_(total_hydrated_volume_);
//    ParallelepipedSeacher_(total_molecular_volume_);
    std::cout << "Total time of execution " << std::setprecision(5)
        << (time(0) - tStart) << " seconds" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}

void SAScore::PipedLoad(std::string file_name) {
    std::ifstream file(file_name);
    Parallelepiped_ temp;
    while (!file.eof()) {
        for (int i = 0; i < 3; ++i)
            file >> temp.coordinates[i];
        file >> temp.length;
        parallelepipeds_.push_back(temp);
    }
    parallelepipeds_.pop_back();
}

void SAScore::CalcIntensity(double Min_q, double Max_q, int N_q, std::complex <double> rho0) {
    std::time_t tStart = time(0);
    rho_0_ = rho0;
    std::cout << "Intensity calculating" << std::endl;;
    IntensityCurveCalc_(Min_q, Max_q, N_q);
    std::cout<< "   Done!!!" << std::endl;
    std::cout << "time of execution "
        << (time(0) - tStart) << " seconds" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}

void SAScore::CalcIntensityBigPiped(double Min_q, double Max_q, int N_q, std::complex <double> rho0) {
    std::time_t tStart = time(0);
    rho_0_ = rho0;
    std::cout << "Intensity calculating" << std::endl;;
    IntensityCurveCalcBigPiped_(Min_q, Max_q, N_q);
    std::cout<< "   Done!!!" << std::endl;
    std::cout << "time of execution "
        << (time(0) - tStart) << " seconds" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
}


void SAScore::CalculatingQ0_() {
    const double pi = acos(-1);

    q_space_.q0_x.resize(q_space_.N_z);
    q_space_.q0_y.resize(q_space_.N_z);
    q_space_.q0_z.resize(q_space_.N_z);

    std::vector <double> theta(q_space_.N_z);
    std::vector <double> phi(q_space_.N_phi);

    for (int i = 0; i < q_space_.N_z; ++i) {
        q_space_.q0_x[i].resize(q_space_.N_phi);
        q_space_.q0_y[i].resize(q_space_.N_phi);
        q_space_.q0_z[i].resize(q_space_.N_phi);
    }

    for (int i = 0; i < q_space_.N_z; ++i) {
        for (int j = 0; j < q_space_.N_phi; ++j) {
            q_space_.q0_z[i][j] = -1.0 + (2.0 / (q_space_.N_z - 1)*i);
        }
        theta[i] = std::acos(q_space_.q0_z[i][0]);
    }

    for (int i = 0; i < q_space_.N_phi; ++i) {
        phi[i] = 2 * pi / q_space_.N_phi*i;
    }

    for (int i = 0; i < q_space_.N_z; ++i) {
        for (int j = 0; j < q_space_.N_phi; ++j) {
            q_space_.q0_x[i][j] = sin(theta[i])*cos(phi[j]);
            q_space_.q0_y[i][j] = sin(theta[i])*sin(phi[j]);
        }
    }
}


void SAScore::IntensityCurveCalc_(double Min_q, double Max_q, int N_q) {
    CalculatingQ0_();
    std::vector <double> q(N_q);
    for (int i = 0; i < N_q; i++) {
        q[i] = Min_q + (Max_q - Min_q) / (N_q - 1)*i;
    }
    int  k = 0;
    for (int i = 0; i < int(N_q / knum_threads) + 1; ++i) {
        bool *thread_state = new bool[knum_threads] {false};
        for (unsigned int j = 0; j < knum_threads; ++j) {
            k = i*knum_threads + j;
            if (k > N_q - 1)
                break;
            thread_[j] = std::thread([this](double q) {
                mtx_.lock();
                std::cout << "|q| = " << q << " Start!!!" << std::endl;
                mtx_.unlock();
                double a = step_ / 2.0;
                std::complex <double> f2 = 0, phi = 0;
                std::complex <double> c = 0, f1 = 0, f = 0;
                double t = 0;
                for (int m = 0; m < q_space_.N_z; ++m) {
                    for (int n = 0; n < q_space_.N_phi; ++n) {
                        f1 = 8 * pow(a, 3)*boost::math::sinc_pi(q * q_space_.q0_x[m][n] * a)*
                            boost_sinc_pi(q * q_space_.q0_y[m][n] * a);

                        f2 = 0;
                        for (auto& piped : parallelepipeds_) {
                            f2 += piped.length * boost::math::sinc_pi(piped.length * a * q *
                                q_space_.q0_z[m][n]) * std::exp(1i* q*
                                (q_space_.q0_x[m][n] * piped.coordinates[int(CubeComponents_::kX)] +
                                    q_space_.q0_y[m][n] * piped.coordinates[int(CubeComponents_::kY)] +
                                    q_space_.q0_z[m][n] * piped.coordinates[int(CubeComponents_::kZ)]));
                        }

                        phi = 0.0;
                        for (auto& iter : pdb_data_) {
                            phi += CalcFactor_(q, iter) * exp(1i*q*
                                (q_space_.q0_x[m][n] * iter[int(PDBComponents_::kX)] +
                                    q_space_.q0_y[m][n] * iter[int(PDBComponents_::kY)] +
                                    q_space_.q0_z[m][n] * iter[int(PDBComponents_::kZ)]));
                        }
                        t += 1;
                        c = phi - rho_0_ * f1 * f2;
                        f += c*std::conj(c);
                    }
                }
                I_.push_back({q, f.real() / t});
                mtx_.lock();
                std::cout << "|q| = " << q << " done!!!" << std::endl;
                mtx_.unlock();
            }, q[k]);
            thread_state[j] = true;
        }

        for (unsigned int j = 0; j < knum_threads; ++j) {
            if (thread_state[j] == true)
                thread_[j].join();
        }
    }
    sort(I_.begin(), I_.end(), ByQ_);

    WriteInFile_(I_, "curve");
}

void SAScore::IntensityCurveCalc_(const std::vector<double>& q) {
    CalculatingQ0_();
    std::vector <double> F;
    size_t N_q = q.size();
    unsigned int  k = 0;
    for (int i = 0; i < int(N_q / knum_threads) + 1; ++i) {
        bool *thread_state = new bool[knum_threads] {false};
        for (unsigned int j = 0; j < knum_threads; ++j) {
            k = i*knum_threads + j;
            if (k > N_q - 1)
                break;
            thread_[j] = std::thread([this](double q) {
                mtx_.lock();
                std::cout << "|q| = " << q << " Start!!!" << std::endl;
                mtx_.unlock();
                double a = step_ / 2.0;
                std::complex <double> f2 = 0, phi = 0;
                std::complex <double> c = 0, f1 = 0, f = 0;
				double t = 0;
                for (int m = 0; m < q_space_.N_z; ++m) {
                    for (int n = 0; n < q_space_.N_phi; ++n) {
                        f1 = 8 * pow(a, 3)*boost::math::sinc_pi(q * q_space_.q0_x[m][n] * a)*
							boost_sinc_pi(q * q_space_.q0_y[m][n] * a);

                        f2 = 0;
                        for (auto& piped : parallelepipeds_) {
                            f2 += piped.length * boost::math::sinc_pi(piped.length * a * q *
                                q_space_.q0_z[m][n]) * std::exp(1i * q *
                                (q_space_.q0_x[m][n] * piped.coordinates[int(CubeComponents_::kX)] +
                                    q_space_.q0_y[m][n] * piped.coordinates[int(CubeComponents_::kY)] +
                                    q_space_.q0_z[m][n] * piped.coordinates[int(CubeComponents_::kZ)]));
                        }

                        phi = 0.0;
                        for (auto& iter : pdb_data_) {
                            phi += CalcFactor_(q, iter)* exp(1i * q *
                                (q_space_.q0_x[m][n] * iter[int(PDBComponents_::kX)] +
                                    q_space_.q0_y[m][n] * iter[int(PDBComponents_::kY)] +
                                    q_space_.q0_z[m][n] * iter[int(PDBComponents_::kZ)]));
                        }
                        t += 1;
                        c = phi - rho_0_ * f1 * f2;
                        f += c*std::conj(c);
                    }
                }
                I_.push_back({q, f.real() / t});
                mtx_.lock();
                std::cout << "|q| = " << q << " done!!!" << std::endl;
                mtx_.unlock();
            }, q[k]);
            thread_state[j] = true;
        }

        for (unsigned int j = 0; j < knum_threads; ++j) {
            if (thread_state[j] == true)
                thread_[j].join();
        }
    }
    sort(I_.begin(), I_.end(), ByQ_);

	WriteInFile_(I_, "curve");
}


double SAScore::GetChi(){
    return chi_;
}

void SAScore::CalculatingChi_(std::vector< std::vector<double> >& data){
        size_t N = data.size();
        std::valarray<double> q(N), Ic(N), sigma(N), Ie(N);

        for (unsigned i = 0; i < N; ++i) {
            q[i] = data[i][0];
            Ie[i] = data[i][1];
            sigma[i] = data[i][2];
            Ic[i] = I_[i][1];
        }
        double c;
        c = (Ie*Ic / (sigma*sigma)).sum() / (Ic*Ic / (sigma*sigma)).sum();
        chi_ = std::sqrt((((Ie - c * Ic) / sigma)*((Ie - c * Ic) / sigma)).sum()/(N-1));

}

void SAScore::CalculatingIntChi(std::vector< std::vector<double> >& data, std::complex <double> rho0){
        std::vector<double> q;
        for (auto i: data){
            q.push_back(i[0]);
        }

        std::time_t tStart = time(0);
        rho_0_ = rho0;
        std::cout << "Intensity calculating" << std::endl;;
        IntensityCurveCalc_(q);
        CalculatingChi_(data);
        std::cout<< "   Done!!!" << std::endl;
        std::cout << "time of execution "
            << (time(0) - tStart) << " seconds" << std::endl;
       std::cout << "------------------------------------------" << std::endl;

}

int SAScore::NumBorderCubes_(std::vector< std::vector<int> >& data){
    unsigned int num = 0;
    for (auto i: data){
        if (i[int(CubeComponents_::kState)]==2)
            ++num;
    }
    return num;
}

double SAScore::ProteinVolume() {
    double volume = 0.0;
    for (auto i:parallelepipeds_)
        volume += i.length;
    return volume*pow(step_, 3);
}

double SAScore::CalcFactor_(double& q, std::vector<double>& atom){
    double f = 0;
    for (unsigned int i = int(PDBComponents_::ka0 ); i <= int(PDBComponents_::ka5); ++i){
        f += atom[i]*std::pow(q, 2*(i-int(PDBComponents_::ka0)));
    }
    return f;
}

double SAScore::RadiusGuinier(){
  const int n1 = 0, n2 = 2;
  return std::sqrt(std::log(I_[n2][1]/I_[n1][1])/(std::pow(I_[n2][0], 2)-std::pow(I_[n1][0],2))*(-3));
}


void SAScore::SetQSpace(int N_z, int N_phi){
    q_space_.N_z = N_z;
    q_space_.N_phi = N_phi;
}


void SAScore::BigPipedParam(){
    double x_min = 0, y_min = 0, z_min = 0,
           x_max = 0, y_max = 0, z_max = 0;

//    for (auto i: pdb_data_){
//        if (x_min > i[int(PDBComponents_::kX)])
//            x_min = i[int(PDBComponents_::kX)];
//        if (x_max < i[int(PDBComponents_::kX)])
//            x_max = i[int(PDBComponents_::kX)];

//        if (y_min > i[int(PDBComponents_::kY)])
//            y_min = i[int(PDBComponents_::kY)];
//        if (y_max < i[int(PDBComponents_::kY)])
//            y_max = i[int(PDBComponents_::kY)];

//        if (z_min > i[int(PDBComponents_::kZ)])
//            z_min = i[int(PDBComponents_::kZ)];
//        if (z_max < i[int(PDBComponents_::kZ)])
//            z_max = i[int(PDBComponents_::kZ)];
//    }


    x_min = big_piped.x[0];
    x_max = big_piped.x[1];

    y_min = big_piped.y[0];
    y_max = big_piped.y[1];

    z_min = big_piped.z[0];
    z_max = big_piped.z[1];

    big_piped_.x_length = x_max - x_min;
    big_piped_.y_length = y_max - y_min;
    big_piped_.z_length = z_max - z_min;

    big_piped_.coordinates[int(CubeComponents_::kX)] = (x_max + x_min) * 0.5;
    big_piped_.coordinates[int(CubeComponents_::kY)] = (y_max + y_min) * 0.5;
    big_piped_.coordinates[int(CubeComponents_::kZ)] = (z_max + z_min) * 0.5;

    big_piped_.volume = big_piped_.x_length * big_piped_.y_length * big_piped_.z_length;
}


//void SAScore::IntensityCurveCalcBigPiped_(double Min_q, double Max_q, int N_q) {
//    BigPipedParam();
//    CalculatingQ0_();
//    std::vector <double> q(N_q);
//    for (int i = 0; i < N_q; i++) {
//        q[i] = Min_q + (Max_q - Min_q) / (N_q - 1)*i;
//    }
//    int  k = 0;
//    for (int i = 0; i < int(N_q / knum_threads) + 1; ++i) {
//        bool *thread_state = new bool[knum_threads] {false};
//        for (unsigned int j = 0; j < knum_threads; ++j) {
//            k = i*knum_threads + j;
//            if (k > N_q - 1)
//                break;
//            thread_[j] = std::thread([this](double q) {
////                mtx_.lock();
////                std::cout << "|q| = " << q << " Start!!!" << std::endl;
////                mtx_.unlock();
//                std::complex <double> c = 0, a = 0, f = 0, phi = 0;
//                double t = 0;
//                for (int m = 0; m < q_space_.N_z; ++m) {
//                   for (int n = 0; n < q_space_.N_phi; ++n){

//                       phi = 0.0;

//                       a = 8 * big_piped_.volume *
//                               boost::math::sinc_pi(q * q_space_.q0_x[m][n] * big_piped_.x_length * 0.5) *
//                               boost::math::sinc_pi(q * q_space_.q0_y[m][n] * big_piped_.y_length * 0.5) *
//                               boost::math::sinc_pi(q * q_space_.q0_z[m][n] * big_piped_.z_length * 0.5) *
//                       std::exp(1i* q*
//                                (q_space_.q0_x[m][n] * big_piped_.coordinates[int(CubeComponents_::kX)] +
//                                 q_space_.q0_y[m][n] * big_piped_.coordinates[int(CubeComponents_::kY)] +
//                                 q_space_.q0_z[m][n] * big_piped_.coordinates[int(CubeComponents_::kZ)]));

//                        for (auto& iter : pdb_data_) {
//                            phi += CalcFactor_(q, iter) * exp(1i*q*
//                                (q_space_.q0_x[m][n] * iter[int(PDBComponents_::kX)] +
//                                    q_space_.q0_y[m][n] * iter[int(PDBComponents_::kY)] +
//                                    q_space_.q0_z[m][n] * iter[int(PDBComponents_::kZ)]));
//                        }
//                        t += 1;
//                        c = phi - rho_0_ * a;
//                        f += c*std::conj(c);
//                    }
//                }
//                I_.push_back({q, f.real() / t});
////                mtx_.lock();
////                std::cout << "|q| = " << q << " done!!!" << std::endl;
////                mtx_.unlock();
//            }, q[k]);
//            thread_state[j] = true;
//        }

//        for (unsigned int j = 0; j < knum_threads; ++j) {
//            if (thread_state[j] == true)
//                thread_[j].join();
//        }
//    }
//    sort(I_.begin(), I_.end(), ByQ_);

//    WriteInFile_(I_, "curve");
//}


void SAScore::IntensityCurveCalcBigPiped_(double Min_q, double Max_q, int N_q) {
    BigPipedParam();
    CalculatingQ0_();
    std::vector <double> Q(N_q);
    for (int i = 0; i < N_q; i++) {
        Q[i] = Min_q + (Max_q - Min_q) / (N_q - 1)*i;
    }

    for(auto q: Q){
                std::cout << "|q| = " << q << " Start!!!" << std::endl;
                std::complex <double> c = 0, a = 0, f = 0, phi = 0;
                double t = 0;
                for (int m = 0; m < q_space_.N_z; ++m) {
                   for (int n = 0; n < q_space_.N_phi; ++n){

                       phi = 0.0;

                       a = 8 * big_piped_.volume *
                               boost::math::sinc_pi(q * q_space_.q0_x[m][n] * big_piped_.x_length * 0.5) *
                               boost::math::sinc_pi(q * q_space_.q0_y[m][n] * big_piped_.y_length * 0.5) *
                               boost::math::sinc_pi(q * q_space_.q0_z[m][n] * big_piped_.z_length * 0.5) *
                       std::exp(1i* q*
                                (q_space_.q0_x[m][n] * big_piped_.coordinates[int(CubeComponents_::kX)] +
                                 q_space_.q0_y[m][n] * big_piped_.coordinates[int(CubeComponents_::kY)] +
                                 q_space_.q0_z[m][n] * big_piped_.coordinates[int(CubeComponents_::kZ)]));

                        for (auto& iter : pdb_data_) {
                            phi += CalcFactor_(q, iter) * exp(1i*q*
                                (q_space_.q0_x[m][n] * iter[int(PDBComponents_::kX)] +
                                    q_space_.q0_y[m][n] * iter[int(PDBComponents_::kY)] +
                                    q_space_.q0_z[m][n] * iter[int(PDBComponents_::kZ)]));
                        }
                        t += 1;
                        c = phi - rho_0_ * a;
                        f += c*std::conj(c);
                    }
                }
                I_.push_back({q, f.real() / t});
                std::cout << "|q| = " << q << " done!!!" << std::endl;
    }

    WriteInFile_(I_, "curve");
}



//void SAScore::IntensityCurveCalcBigPiped_(double Min_q, double Max_q, int N_q) {
//    BigPipedParam();
//    CalculatingQ0_();
//    std::vector <double> Q(N_q);
//    std::vector <std::thread> threads;
//    for (int i = 0; i < N_q; i++) {
//        Q[i] = Min_q + (Max_q - Min_q) / (N_q - 1)*i;
//    }

//    for(double &q: Q){
//        threads.push_back( std::thread([this](double q){
//                std::cout << "|q| = " << q << " Start!!!" << std::endl;
//                std::complex <double> c = 0, a = 0, f = 0, phi = 0;
//                double t = 0;
//                for (int m = 0; m < q_space_.N_z; ++m) {
//                   for (int n = 0; n < q_space_.N_phi; ++n){

//                       phi = 0.0;

//                       a = 8 * big_piped_.volume *
//                               boost::math::sinc_pi(q * q_space_.q0_x[m][n] * big_piped_.x_length * 0.5) *
//                               boost::math::sinc_pi(q * q_space_.q0_y[m][n] * big_piped_.y_length * 0.5) *
//                               boost::math::sinc_pi(q * q_space_.q0_z[m][n] * big_piped_.z_length * 0.5) *
//                       std::exp(1i* q*
//                                (q_space_.q0_x[m][n] * big_piped_.coordinates[int(CubeComponents_::kX)] +
//                                 q_space_.q0_y[m][n] * big_piped_.coordinates[int(CubeComponents_::kY)] +
//                                 q_space_.q0_z[m][n] * big_piped_.coordinates[int(CubeComponents_::kZ)]));

//                        for (auto& iter : pdb_data_) {
//                            phi += CalcFactor_(q, iter) * exp(1i*q*
//                                (q_space_.q0_x[m][n] * iter[int(PDBComponents_::kX)] +
//                                    q_space_.q0_y[m][n] * iter[int(PDBComponents_::kY)] +
//                                    q_space_.q0_z[m][n] * iter[int(PDBComponents_::kZ)]));
//                        }
//                        t += 1;
//                        c = phi - rho_0_ * a;
//                        f += c*std::conj(c);
//                    }
//                }
//                I_.push_back({q, f.real() / t});
//                std::cout << "|q| = " << q << " done!!!" << std::endl;
//        },q));
//    }

//    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

//    sort(I_.begin(), I_.end(), ByQ_);

//    WriteInFile_(I_, "curve");
//}
