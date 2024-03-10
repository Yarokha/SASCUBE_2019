#include "PDBParser.h"


int PDBParser::ConverterNum_(std::string str) {
    int const_shift = 0;
    int num = 100000;
    if ((str[0] - 55) > 9) {
        num += int(str[0] - 65)*std::pow(36, 4);

        for (unsigned int i = 1; i < 5; ++i) {
            if (str[i] - 48 > 9) const_shift = 55;
            else const_shift = 48;
            num += (str[i] - const_shift) *std::pow(36, 4 - i);

        }
    }
    else {
        std::stringstream temp_io;
        temp_io << str;
        temp_io >> num;
        temp_io.clear();
    }
    return num;
}

PDBParser::ATOM PDBParser::AtomParse_(std::string atom_string) {

    ATOM atom_temp;
    std::stringstream temp_io;

    temp_io << atom_string.substr(0, 6);
    temp_io >> atom_temp.type;
    temp_io.clear();

    atom_temp.serial = ConverterNum_(atom_string.substr(6, 5));

    temp_io << atom_string.substr(12, 4);
    temp_io >> atom_temp.name;
    temp_io.clear();

    atom_temp.altLoc = atom_string[16];

    temp_io << atom_string.substr(17, 3);
    temp_io >> atom_temp.resName;
    temp_io.clear();

    atom_temp.chainID = atom_string[21];

    atom_temp.resSeq = atom_string.substr(22, 4);

    atom_temp.iCode = atom_string[26];

    temp_io << atom_string.substr(30, 8);
    temp_io >> atom_temp.x;
    temp_io.clear();

    temp_io << atom_string.substr(38, 8);
    temp_io >> atom_temp.y;
    temp_io.clear();

    temp_io << atom_string.substr(46, 8);
    temp_io >> atom_temp.z;
    temp_io.clear();

    temp_io << atom_string.substr(54, 6);
    temp_io >> atom_temp.occupancy;
    temp_io.clear();

    temp_io << atom_string.substr(60, 6);
    temp_io >> atom_temp.tempFactor;
    temp_io.clear();

    temp_io << atom_string.substr(76, 2);
    temp_io >> atom_temp.element;
    temp_io.clear();

    temp_io << atom_string.substr(78, 2);
    temp_io >> atom_temp.charge;
    temp_io.clear();

    return atom_temp;
}


void PDBParser::ParsePDB_() {
	std::ifstream pdb_file(path_to_pdb_);
	std::string atom_line;
	ATOM atom_temp;
	while (std::getline(pdb_file, atom_line)) {
		std::stringstream temp_io;
		std::string temp_input, check, amino_temp;

		for (unsigned int i = 0; i < 6; i++)
			temp_input += atom_line[i];
		temp_io << temp_input;
		temp_io >> check;
		temp_io.clear();
		if (check == "HEADER" && PDB_id_ == "") {
			PDB_id_ = atom_line.substr(62, 4);
		}
		else if (check == "TITLE" && protein_title_ == "") {
			protein_title_ = atom_line.substr(10, atom_line.size());
			std::clog << protein_title_ << std::endl
				<< "PDB ID: " << PDB_id_ << std::endl;
		}
		else if (check == "ATOM" || check == "HETATM") {
			atom_temp = AtomParse_(atom_line);
			amino_temp = aminoacid_[atom_temp.resName][atom_temp.name];
			if (amino_temp == "") {
				std::clog << "No information about" << std::endl
					<< "Atom type: " << check << std::endl
					<< "Residue Name: " << atom_temp.resName << std::endl
					<< "Atom name: " << atom_temp.name << std::endl
					<< "--------------------" << std::endl;
			}
			atom_temp.radius = &aminoacid_info_[amino_temp].radius;
			atom_temp.formFactor = &aminoacid_info_[amino_temp].formFactor;
			atom_temp.formFactorCoef = &aminoacid_info_[amino_temp].formFactorCoef;
			if (check == "ATOM") {
				++chains_[atom_temp.chainID].number_atom;
			}
			else {
				if (atom_temp.resName == "HOH")
					++chains_[atom_temp.chainID].number_hetatom_hoh;
				else
					++chains_[atom_temp.chainID].number_hetatom_nonhoh;
			}
			protein_data_.push_back(atom_temp);
			
		}


        //experimental
        else if (check == "REMARK") {

            //REMARK 250 MIN_X = -27.315 MAX_X = 23.047 DIFFERENCE_X = 50.362
            //REMARK 250 MIN_Y = -27.048 MAX_Y = 24.886 DIFFERENCE_Y = 51.934
            //REMARK 250 MIN_Z = -23.432 MAX_Z = 27.084 DIFFERENCE_Z = 50.516

            if (atom_line.substr(7, 9) == "250 MIN_X") {
                BigParallelepiped_(atom_line, big_piped.x);
            }
            else if (atom_line.substr(7, 9) == "250 MIN_Y") {
                BigParallelepiped_(atom_line, big_piped.y);
            }
            else if (atom_line.substr(7, 9) == "250 MIN_Z") {
                BigParallelepiped_(atom_line, big_piped.z);
        }
        }

	}
	pdb_file.close();
}


void PDBParser::Aminocid_() {
	std::ifstream amino_file(path_to_aminoacid_inf_);
    std::vector<std::string> temp;
    std::string temp_line, temp_amino;
    double coef = 0;
	
	while (std::getline(amino_file, temp_line)){
		std::istringstream line(temp_line);
        line >> temp_amino;
        line >> aminoacid_info_[temp_amino].formFactor;
        line >> aminoacid_info_[temp_amino].radius;
        while (line){
            line>>coef;
            aminoacid_info_[temp_amino].formFactorCoef.push_back(coef);
        }
        aminoacid_info_[temp_amino].formFactorCoef.pop_back();
	}
	amino_file.close();

    std::string filename;
    QDir directory(QString::fromStdString(path_to_aminoacid_folder_));
    QStringList AminoFiles = directory.entryList(QStringList() << "*.txt" << "*.txt",QDir::Files);
    foreach(QString Filename, AminoFiles) {

        filename = Filename.toStdString();
        filename = filename.substr(0, filename.size() - 4);
        std::ifstream amino_file(path_to_aminoacid_folder_+'/'+Filename.toStdString());
        std::string amino_line, amino;

        while (std::getline(amino_file, amino_line)) {
            std::istringstream line(amino_line);

            while (line) {
                line >> amino;
                temp.push_back(amino);
            }

            aminoacid_[filename][temp[0]] = temp[1];
            temp.clear();
        }
    }

}


unsigned int PDBParser::AtomCount()const {
	unsigned int count = 0;
	for (auto i : chains_)
		count += (i.second.number_atom + i.second.number_hetatom_nonhoh);
	return count;
}


void PDBParser::BigParallelepiped_(std::string line, std::vector<double> &a) {
    std::stringstream temp_io;

    temp_io << line.substr(17, 7);
    temp_io >> a[0];
    a[0] *= 0.1;
    temp_io.clear();

    temp_io << line.substr(31, 7);
    temp_io >> a[1];
    a[1] *= 0.1;
    temp_io.clear();
    std::cout << line << std::endl;

}


PDBParser::PDBParser(std::string path_to_pdb, std::string path_to_aminoacid_folder, std::string path_to_aminoacid_inf):
    path_to_pdb_(path_to_pdb), path_to_aminoacid_folder_(path_to_aminoacid_folder), path_to_aminoacid_inf_(path_to_aminoacid_inf){
    Aminocid_();
	ParsePDB_();
    std::clog << "Detected " << ChainsNumber() << " chain(s)" << std::endl <<
		"Atom Count: " << AtomCount() << std::endl;
}


PDBParser::~PDBParser(){
}
