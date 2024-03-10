#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <QDir>

class PDBParser
{
private:
	std::string protein_title_, PDB_id_;

	struct ATOM {
		//	COLUMNS		DATA			TYPE		FIELD DEFINITION
		//	------------------------------------------------------------------------------------ -
		//	1 - 6		Record			name		"ATOM "
		//	7 - 11		Integer			serial		Atom serial number.
		//	13 - 16		Atom			name		Atom name.
		//	17			Character		altLoc		Alternate location indicator.
		//	18 - 20		Residue	name	resName		Residue name.
		//	22			Character		chainID		Chain identifier.
		//	23 - 26		Integer			resSeq		Residue sequence number.
		//	27			AChar			iCode		Code for insertion of residues.
		//	31 - 38		Real(8.3)		x			Orthogonal coordinates for X in Angstroms.
		//	39 - 46		Real(8.3)		y			Orthogonal coordinates for Y in Angstroms.
		//	47 - 54		Real(8.3)		z			Orthogonal coordinates for Z in Angstroms.
		//	55 - 60		Real(6.2)		occupancy	Occupancy.
		//	61 - 66		Real(6.2)		tempFactor	Temperature factor.
		//	77 - 78		LString(2)		element		Element symbol, right - justified.
		//	79 - 80		LString(2)		charge		Charge on the atom.

		std::string type;
        int serial;
		std::string name;
		char altLoc;
		std::string resName;
        char chainID;
        std::string resSeq;
		char iCode;
		double x, y, z, occupancy, tempFactor;
		std::string element, charge;
		//-------------------------------------------------------------------------
		//additional information
		double *radius;
		double *formFactor;
		//-------------------------------------------------------------------------
		//formFactor coefficients
		std::vector<double> *formFactorCoef;
	};
	struct AMINO {
		double radius;
		double formFactor;
		std::vector<double> formFactorCoef;
	};
	struct CHAIN {
		unsigned int number_atom = 0, number_hetatom_hoh = 0 , number_hetatom_nonhoh = 0;
    };

	std::map < std::string, AMINO > aminoacid_info_;
	std::map < std::string, std::map < std::string, std::string> > aminoacid_;
	std::map < char, CHAIN > chains_;
	std::vector<ATOM> protein_data_;

	std::string path_to_pdb_, path_to_aminoacid_inf_, path_to_aminoacid_folder_;
	
	ATOM AtomParse_(std::string atom_string);
	void ParsePDB_();
	void Aminocid_();
    int ConverterNum_(std::string str);
    void BigParallelepiped_(std::string line, std::vector<double> &a);

public:
	PDBParser(std::string path_to_pdb, std::string path_to_aminoacid_folder, std::string path_to_aminoacid_inf);
	~PDBParser();
	const std::vector<ATOM>& ProteinData()const { return protein_data_; }
	const std::string& ProteinTitle()const { return protein_title_; }
	const std::string& PDBId()const { return PDB_id_; }
	const std::map < char, CHAIN >& ChainInfo()const { return chains_; }
    size_t ChainsNumber() { return chains_.size(); }
    unsigned int AtomCount() const;
    struct PARALLELEPIPED {
            std::vector<double>
                x = std::vector<double>(2),
                y = std::vector<double>(2),
                z = std::vector<double>(2);
    } big_piped;
};

