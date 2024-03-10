#include "SAXSParser.h"


std::vector < std::vector<double> > SAXSParser::AtomTable(bool serial, bool nm, bool hetatom, bool hoh_hetatom, bool hem_hetatom) const {
	std::vector< std::vector<double> > atom_table;
	
    double k = 0.1;
	if (!nm)
		k = 1.0;
	
    for (auto& i : ProteinData()) {
        if ((i.resName == "HOH" && hoh_hetatom == false)|| (i.resName == "HEM" && hem_hetatom == false) || (i.type == "HETATM" && i.resName != "HEM" && i.resName != "HOH"
                                                                                                            && hetatom == false))
			continue;
		if (serial)
			atom_table.push_back({double(i.serial), i.x*k, i.y*k, i.z*k, *i.radius*k, *i.formFactor});
		else
			atom_table.push_back({i.x*k, i.y*k, i.z*k, *i.radius*k, *i.formFactor});
		
		for (auto& coef : *i.formFactorCoef)
			atom_table.back().push_back(coef);

        if ((atom_table.back().size() < 12 && serial) || (atom_table.back().size() < 11 && !serial))
                    atom_table.pop_back();
	}

	return atom_table;
}

void SAXSParser::SaveDataToFile(bool serial, bool nm, bool hetatom, bool hoh_hetatom, bool hem_hetatom, std::string file_path) const{
    std::ofstream datafile(file_path + '/' + PDBId()+"_data.txt");
    for (auto& atom : AtomTable(serial, nm, hetatom, hoh_hetatom, hem_hetatom)) {
		for (auto& i : atom) {
			datafile << i << '\t';
		}
		datafile << std::endl;
	}
	datafile.close();
}

SAXSParser::SAXSParser(std::string path_to_pdb, std::string path_to_aminoacid_folder, std::string path_to_aminoacid_inf):
	PDBParser(path_to_pdb, path_to_aminoacid_folder, path_to_aminoacid_inf)
{
}


SAXSParser::~SAXSParser()
{
}
