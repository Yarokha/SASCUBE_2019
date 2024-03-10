#pragma once
#include "PDBParser.h"

class SAXSParser :
	public PDBParser
{
public:
	std::vector < std::vector<double> > AtomTable(bool serial = false, bool nm = true/*true in nanometer, false in Angstrom */,
        bool hetatom = false, bool hoh_hetatom = false, bool hem_hetatom = false) const;
    void SaveDataToFile(bool serial, bool nm, bool hetatom, bool hoh_hetatom, bool hem_hetatom, std::string file_path) const;
	SAXSParser(std::string path_to_pdb, std::string path_to_aminoacid_folder, std::string path_to_aminoacid_inf);
	~SAXSParser();
};

