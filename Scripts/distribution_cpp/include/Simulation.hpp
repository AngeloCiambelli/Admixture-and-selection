#ifndef DEF_SIMULATION_HPP
#define DEF_SIMULATION_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <string>
#include <Python.h>

#include "matplotlibcpp.h" // Display with matplotlib
#include "Vecteur.hpp"
#include "Admixture.hpp"

namespace plt = matplotlibcpp;
using namespace std;

//===========================================================================
//                          Description
//===========================================================================
//
// Class Simulation : Simulate addmixture with random mating and if wanted selection 
// given a contribution vector (contribution) and a selection matrix (Phi).
// 
//
//===========================================================================
//                          Define Simulation class
//===========================================================================

class Simulation
{
    public:
    Admixture final_admixture_proportion;  
    Vecteur<Vecteur<float>> contribution; //List of contribution
    Vecteur<Vecteur<Vecteur<float>>> Phi; //List of selection matrices
    int nPop; 
    int nGen;
    string fileName;

    //n source pop
    Simulation(){}
    Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int n, int gen, string filename, bool allplot, bool allsave, bool endplot, bool endsave); //Simulation avec panmixie+selection
    Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int n, int gen, string filename);

    //Two source pop
    Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int gen, string filename, bool allplot, bool allsave, bool endplot, bool endsave); //Simulation avec panmixie+selection
    Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int gen, string filename);

    //Discretized two source pop
    Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int gen, string filename, int nStep, bool allplot, bool allsave, bool endplot, bool endsave); //Simulation avec panmixie+selection
    Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int gen, string filename, int nStep);

    Simulation(const Vecteur<Vecteur<float>> &cont, int gen, string filename); //Simulation avec panmixie seule
    Simulation(const Vecteur<Vecteur<float>> &cont, int gen, string filename, int nStep);
};


//===========================================================================
//                          Constructors
//===========================================================================

Simulation::Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int n, int gen, string filename, bool allplot, bool allsave, bool endplot, bool endsave)
{ 
    //Save parameters
    contribution = cont;
    Phi = phi;
    nGen = gen;
    nPop = n;
    fileName = filename;

    //Base vectors
    Vecteur<Vecteur<float>> e(nPop, Vecteur<float> (nPop, 0));
    for (int k=0; k<nPop; k++)
    {
        e[k][k]= 1;
    }

    //Save parameters in a .txt file
    ofstream file(fileName+"parameters.txt");
    file << "Phi:" << Phi << endl;
    file << "Contribution:" << contribution << endl;
    file.close();

    //Produce simulation
    Admixture pool;

    for (int i=0; i<nGen; i++)
    {   
        Admixture pool2;
        Admixture pool3;

        pool.RandomMating(pool2, contribution[i], e);
        pool2.Selection(pool3, Phi[i]);
        pool = pool3;

        cout << "generation : " << i+1 << endl  << pool << endl << endl;

        if (allplot == 1)
        {
            pool.histogram(fileName+"H_1_generation_"+to_string(i+1)+".png", false, 1, pow(2,i+1)+1, i);
        }
        
        if (allsave == 1)
        {
            pool.write_to_txt_R(i+1, fileName, nPop);
        }
    }
    
    final_admixture_proportion = pool;

    if (endplot == 1)
    {
        final_admixture_proportion.histogram(fileName+"H_1_generation_"+to_string(nGen)+".png", false, 1, pow(2,nGen)+1, nGen-1);
    }

    if (endsave == 1)
    {
        final_admixture_proportion.write_to_txt_R(nGen, fileName, nPop);

        //Create .txt for summary statistics of the N simulations referenced by their id
        ofstream newfile;
        newfile.open(fileName+"sumStat.txt");
        newfile << "id;mean;variance;skewness;kurtosis;mode;min;max;Q1;Q2;Q3;Q4;Q5;Q6;Q7;Q8;Q9\n";
        map<string, float> sumStat = final_admixture_proportion.summaryStat(10);
        string id = "to_find_sim";
        newfile << id << ";" << sumStat["mean"] << ";" 
                << sumStat["variance"] << ";" 
                << sumStat["skewness"] << ";" 
                << sumStat["kurtosis"] << ";" 
                << sumStat["mode"] << ";" 
                << sumStat["min"] << ";" 
                << sumStat["max"] << ";" 
                << sumStat["Q1"] << ";" 
                << sumStat["Q2"] << ";" 
                << sumStat["Q3"] << ";" 
                << sumStat["Q4"] << ";" << sumStat["Q5"] 
                << ";" << sumStat["Q6"] << ";" 
                << sumStat["Q7"] << ";" 
                << sumStat["Q8"] << ";" 
                << sumStat["Q9"] << "\n";
        newfile.close();
    }
}

Simulation::Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int n, int gen, string filename)
{ 
    //Save parameters
    contribution = cont;
    Phi = phi;
    nGen = gen;
    nPop = n;
    fileName = filename;

    //Base vectors
    Vecteur<Vecteur<float>> e(nPop, Vecteur<float> (nPop, 0));
    for (int k=0; k<nPop; k++)
    {
        e[k][k]= 1;
    }

    //Save parameters in a .txt file
    ofstream file(fileName+"parameters.txt");
    file << "Phi:" << Phi << endl;
    file << "Contribution:" << contribution << endl;
    file.close();

    //Produce simulation
    Admixture pool;

    for (int i=0; i<nGen; i++)
    {   
        Admixture pool2;
        Admixture pool3;

        pool.RandomMating(pool2, contribution[i], e);
        pool2.Selection(pool3, Phi[i]);
        pool = pool3;
        
        //cout << "generation : " << i+1 << endl ; //<< pool << endl << endl;
    }
    
    final_admixture_proportion = pool;
}


Simulation::Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int gen, string filename, bool allplot, bool allsave, bool endplot, bool endsave)
{ 
    //Save parameters
    contribution = cont;
    Phi = phi;
    nGen = gen;
    nPop = 2;
    fileName = filename;

    //Save parameters in a .txt file
    ofstream file(fileName+"parameters.txt");
    file << "Phi:" << Phi << endl;
    file << "Contribution:" << contribution << endl;
    file.close();

    //Produce simulation
    Admixture pool;

    for (int i=0; i<nGen; i++)
    {  
        Admixture pool2;
        Admixture pool3;

        pool.RandomMating_H1(pool2, contribution[i]);
        pool2.Selection_H1(pool3, Phi[i]);
        pool = pool3;

        //cout << "generation : " << i+1 << endl << pool << endl << endl;
        Admixture converted_pool = pool.convert_fraction_H1_to_fraction();

        if ((allplot == 1) and ((i==0) or (i==3) or (i==7) or (i==11) or (i==15)))
        {   
            cout << i << endl;
            converted_pool.histogram(fileName+"H_1_generation_"+to_string(i+1)+".pdf", false, 1, pow(2,(nGen))+1, i);
        }
        
        if (allsave == 1)
        {
            converted_pool.write_to_txt_R(i+1, fileName, nPop);
        }

    }
    
    final_admixture_proportion = pool;

    Admixture converted_pool = pool.convert_fraction_H1_to_fraction();

    if (endplot == 1)
    {
        converted_pool.histogram(fileName+"H_1_generation_"+to_string(nGen)+".png", false, 1, pow(2,nGen)+1, nGen-1);
    }

    if (endsave == 1)
    {
        converted_pool.write_to_txt_R(nGen, fileName, nPop);

        //Create .txt for summary statistics of the N simulations referenced by their id
        ofstream newfile;
        newfile.open(fileName+"sumStat.txt");
        newfile << "id;mean;variance;skewness;kurtosis;mode;min;max;Q1;Q2;Q3;Q4;Q5;Q6;Q7;Q8;Q9\n";
        map<string, float> sumStat = converted_pool.summaryStat(10);
        string id = "to_find_sim";
        newfile << id << ";" << sumStat["mean"] << ";" 
                << sumStat["variance"] << ";" 
                << sumStat["skewness"] << ";" 
                << sumStat["kurtosis"] << ";" 
                << sumStat["mode"] << ";" 
                << sumStat["min"] << ";" 
                << sumStat["max"] << ";" 
                << sumStat["Q1"] << ";" 
                << sumStat["Q2"] << ";" 
                << sumStat["Q3"] << ";" 
                << sumStat["Q4"] << ";" << sumStat["Q5"] 
                << ";" << sumStat["Q6"] << ";" 
                << sumStat["Q7"] << ";" 
                << sumStat["Q8"] << ";" 
                << sumStat["Q9"] << "\n";
        newfile.close();
    }
}

Simulation::Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int gen, string filename)
{ 
    //Save parameters
    contribution = cont;
    Phi = phi;
    nGen = gen;
    nPop = 2;
    fileName = filename;

    //Save parameters in a .txt file
    ofstream file(fileName+"parameters.txt");
    file << "Phi:" << Phi << endl;
    file << "Contribution:" << contribution << endl;
    file.close();

    //Produce simulation
    Admixture pool;

    for (int i=0; i<nGen; i++)
    {   
        Admixture pool2;
        Admixture pool3;
        
        pool.RandomMating_H1(pool2, contribution[i]);
        pool2.Selection_H1(pool3, Phi[i]);
        
        pool = pool3;
    }
    
    final_admixture_proportion = pool;
}

Simulation::Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int gen, string filename, int nStep, bool allplot, bool allsave, bool endplot, bool endsave)
{ 
    //Save parameters
    contribution = cont;
    Phi = phi;
    nGen = gen;
    nPop = 2;
    fileName = filename;

    Vecteur<float> discretization(nStep,0);
    for (int i=0; i<nStep ;i++)
    {
        discretization[i] = float(i)/(nStep-1);
    }


    //Save parameters in a .txt file
    ofstream file(fileName+"parameters.txt");
    file << "Phi:" << Phi << endl;
    file << "Contribution:" << contribution << endl;
    file.close();

    //Produce simulation
    Admixture pool;

    for (int i=0; i<nGen; i++)
    {  
        Admixture pool2;
        Admixture pool3;

        pool.RandomMating_H1_discrete(pool2, contribution[i], discretization);
        pool2.Selection_H1(pool3, Phi[i]);
        pool = pool3;

        Admixture converted_pool = pool.convert_fraction_H1_to_fraction();
        //cout << "generation : " << i+1 << endl << converted_pool << endl << endl; 

        if (allplot == 1)
        {
            converted_pool.histogram(fileName+"H_1_generation_"+to_string(i+1)+".pdf", false, 1, nStep, i);
        }
        
        if (allsave == 1)
        {
            converted_pool.write_to_txt_R(i+1, fileName, nPop);
        }

    }
    
    final_admixture_proportion = pool;

    Admixture converted_pool = pool.convert_fraction_H1_to_fraction();

    if (endplot == 1)
    {
        converted_pool.histogram(fileName+"H_1_generation_"+to_string(nGen)+".png", false, 1, nStep, nGen-1);
    }

    if (endsave == 1)
    {
        converted_pool.write_to_txt_R(nGen, fileName, nPop);

        //Create .txt for summary statistics of the N simulations referenced by their id
        ofstream newfile;
        newfile.open(fileName+"sumStat.txt");
        newfile << "id;mean;variance;skewness;kurtosis;mode;min;max;Q1;Q2;Q3;Q4;Q5;Q6;Q7;Q8;Q9\n";
        map<string, float> sumStat = converted_pool.summaryStat(10);
        string id = "to_find_sim";
        newfile << id << ";" << sumStat["mean"] << ";" 
                << sumStat["variance"] << ";" 
                << sumStat["skewness"] << ";" 
                << sumStat["kurtosis"] << ";" 
                << sumStat["mode"] << ";" 
                << sumStat["min"] << ";" 
                << sumStat["max"] << ";" 
                << sumStat["Q1"] << ";" 
                << sumStat["Q2"] << ";" 
                << sumStat["Q3"] << ";" 
                << sumStat["Q4"] << ";" << sumStat["Q5"] 
                << ";" << sumStat["Q6"] << ";" 
                << sumStat["Q7"] << ";" 
                << sumStat["Q8"] << ";" 
                << sumStat["Q9"] << "\n";
        newfile.close();
    }
}

Simulation::Simulation(const Vecteur<Vecteur<float>> &cont, const Vecteur<Vecteur<Vecteur<float>>> &phi, int gen, string filename, int nStep)
{ 
    //Save parameters
    contribution = cont;
    Phi = phi;
    nGen = gen;
    nPop = 2;
    fileName = filename;

    Vecteur<float> discretization(nStep,0);
    for (int i=0; i<nStep ;i++)
    {
        discretization[i] = float(i)/(nStep-1);
    }

    //Save parameters in a .txt file
    ofstream file(fileName+"parameters.txt");
    file << "Phi:" << Phi << endl;
    file << "Contribution:" << contribution << endl;
    file.close();

    //Produce simulation
    Admixture pool;

    for (int i=0; i<nGen; i++)
    {   
        Admixture pool2;
        Admixture pool3;
        
        pool.RandomMating_H1_discrete(pool2, contribution[i], discretization);
        pool2.Selection_H1(pool3, Phi[i]);
        
        pool = pool3;
    }
    
    final_admixture_proportion = pool;
}

Simulation::Simulation(const Vecteur<Vecteur<float>> &cont, int gen, string filename)
{ 
    //Save parameters
    contribution = cont;
    nGen = gen;
    nPop = 2;
    fileName = filename;

    //Save parameters in a .txt file
    ofstream file(fileName+"parameters.txt");
    file << "Contribution:" << contribution << endl;
    file.close();

    //Produce simulation
    Admixture pool;

    for (int i=0; i<nGen; i++)
    {   
        Admixture pool2;
        //cout << i << endl;
        pool.RandomMating_H1(pool2, contribution[i]);
        pool=pool2;
    }

    final_admixture_proportion = pool;
}

Simulation::Simulation(const Vecteur<Vecteur<float>> &cont, int gen, string filename, int nStep)
{ 
    //Save parameters
    contribution = cont;
    nGen = gen;
    nPop = 2;
    fileName = filename;

    //MAke discretization point list
    Vecteur<float> discretization(nStep,0);
    for (int i=0; i<nStep ;i++)
    {
        discretization[i] = float(i)/(nStep-1);
    }

    //Save parameters in a .txt file
    ofstream file(fileName+"parameters.txt");
    file << "Contribution:" << contribution << endl;
    file.close();

    //Produce simulation
    Admixture pool;

    for (int i=0; i<nGen; i++)
    {   
        Admixture pool2;
        
        pool.RandomMating_H1_discrete(pool2, contribution[i], discretization);
        
        pool = pool2;
    }
    
    final_admixture_proportion = pool;
}


#endif