#ifndef DEF_ADMIXTURE_HPP
#define DEF_ADMIXTURE_HPP

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "Vecteur.hpp"
#include "matplotlibcpp.h" // Display with matplotlib
#include <Python.h>
#include <map>

namespace plt = matplotlibcpp;
using namespace std;

//===========================================================================
//                          Description
//===========================================================================
//
// Class Admixture
//
//===========================================================================
//                          Define Admixture class
//===========================================================================


class Admixture 
{
    public:
    map<Vecteur<float>, float> fraction;
    map<float, float> fraction_H1;
    Vecteur<float> discretization; //Sorted vector of discretization grid of [0,1]

    Admixture(){};
    Admixture(Vecteur<Vecteur<float>> e, Vecteur<float> p);

    void histogram(string figName, bool show, int nPop, int nStep, int i); 
    void write_to_txt_R(int gen, string filename, int nPop);
    map<string, float> summaryStat(int nQ = 10);
    Vecteur<float> firstMoments(int N);
    Admixture convert_fraction_H1_to_fraction();
    
    Admixture& RandomMating(Admixture &H_t_tilde_plus_1, Vecteur<float> contribution, Vecteur<Vecteur<float>> e);
    Admixture& RandomMating_H1(Admixture &H_t_tilde_plus_1, Vecteur<float> contribution);
    Admixture& RandomMating_H1_discrete(Admixture &H_t_tilde_plus_1, Vecteur<float> contribution, Vecteur<float> discretization);

    Admixture& Selection(Admixture &H_t ,Vecteur<Vecteur<float>> Phi);
    Admixture& Selection_H1(Admixture &H_t ,Vecteur<Vecteur<float>> Phi);
};

//===========================================================================
//                          Member functions
//===========================================================================

Admixture::Admixture(Vecteur<Vecteur<float>> e, Vecteur<float> p)
{
    if (e.size()!=p.size()){cout << "Size btw gen and proba" << endl; exit(1);}
    
    for (int k=0; k<e.size();k++)
    {
        if (p[k]>1){cout << "Probability greater than 1" << endl; exit(1);}
        fraction[e[k]] = p[k];
    }
}

void Admixture::histogram(string figName, bool show, int nPop, int nStep, int i) 
{
    // Prepare the data: x = Admixture_k(t) values and weights = Admixture_k(t) probability
    std::vector<float> x;
    std::vector<float> weights;
 
    for (const auto& value : this->fraction) {
        x.push_back(value.first[0]);
        weights.push_back(value.second);
    }

    // Convert x and weights to Python-compatible strings
    std::ostringstream xStream, weightsStream;
    xStream << "[";
    weightsStream << "[";
    for (size_t i = 0; i < x.size(); ++i) {
        xStream << x[i];
        weightsStream << weights[i];
        if (i < x.size() - 1) {
            xStream << ", ";
            weightsStream << ", ";
        }
    }
    xStream << "]";
    weightsStream << "]";

    // Inject data into the Python script
    std::ostringstream python_code;
    python_code << R"(
import matplotlib.pyplot as plt
import numpy as np
import os

# Ensure the directory exists
os.makedirs(os.path.dirname(')" << figName << R"('), exist_ok=True)

def plot_histogram(data, weights, fig_name, show, nStep, i):

    # Plot histogram using plt.hist with weights and bins
    plt.hist(
        data,
        bins=nStep,
        weights=weights,
        color='blue',
        edgecolor='black',
        align='mid',
    )
    #plt.xlabel("$H_1(t)$ ", fontsize=18)
    #plt.ylabel("$P(H_1(t))$ ", fontsize=18)
    plt.xlim(-0.1, 1.1)
    plt.ylim(0, 0.8/(1.4**i))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    if show:
        plt.show()
    plt.gcf().set_size_inches(8.5, 5)
    plt.savefig(fig_name)
    plt.clf()

# Plot the histogram
plot_histogram()" << xStream.str() << ", " << weightsStream.str() << ", '" << figName << "', " << (show ? "True" : "False") << ", " << nStep << ", " << i << ")";
        
    // Execute the Python script
    PyRun_SimpleString(python_code.str().c_str());
}

void Admixture::write_to_txt_R(int gen, string filename, int nPop)
{
    ofstream myfile;
    myfile.open(filename+"addmixed_proportions_generation_"+to_string(gen)+".txt");

    for (int i=0; i<nPop; i++)
    {
        myfile << "prop_s" + to_string(i+1) + ";";
    }
    myfile << "proba \n";

    for (auto j = (this->fraction).begin(); j != (this->fraction).end(); j++)
    {   
        for (int i=0; i<nPop; i++)
        {
            myfile << j->first[i] << ";";
        }
        myfile << j->second << "\n";
    }
    myfile.close();
}

map<string, float> Admixture::summaryStat(int nQ)
{
    //Final stat of interest
    map<string, float> sumStat;
    sumStat["min"] = 1.f, sumStat["max"] = 0.f, sumStat["mean"] = 0.f, sumStat["variance"] = 0.f, sumStat["skewness"] = 0.f, sumStat["kurtosis"] = 0.f, sumStat["mode"] = 0.f;
    
    // Quantiles (stored in the map with key format: "Q1", "Q2", ...)
    std::vector<std::pair<std::string, float>> quantileKeys;
    for (int i = 1; i < nQ; ++i) {
        quantileKeys.push_back({"Q" + std::to_string(i), i / static_cast<float>(nQ)});
    }

    //Intermidiate variable
    float weightSum = 0.f, weightSumSq = 0.f, weightSumCube = 0.f, weightSumQuart = 0.f, min_loc = 1.f, max_loc = 0.f, max_proba = 0.f;
    vector<pair<float, float>> sortedData;

    for (auto j = (this->fraction).begin(); j != (this->fraction).end(); j++)
    {
        weightSum += j->second * j->first[0];
        min_loc = min(min_loc, j->first[0]);
        max_loc = max(max_loc, j->first[0]);

        if (j->second > max_proba)
        {
            max_proba = j->second;
            sumStat["mode"] = j->first[0];
        }

        //Save the value and proba in vector for easier processing
        sortedData.emplace_back(j->first[0], j->second);
    }

    for (auto j = (this->fraction).begin(); j != (this->fraction).end(); j++)
    {
        float diff = j->first[0] - weightSum;
        weightSumSq += j->second * pow(diff, 2);
        weightSumCube += j->second * pow(diff, 3);
        weightSumQuart += j->second * pow(diff, 4);
    }

    sumStat["max"] = max_loc;
    sumStat["min"] = min_loc; 
    sumStat["mean"] = weightSum;
    sumStat["variance"] = weightSumSq;
    sumStat["skewness"] = weightSumCube / pow(sumStat["variance"], 1.5);
    sumStat["kurtosis"] = weightSumQuart / (sumStat["variance"] * sumStat["variance"]);

    //Sort data for quantiles
    std::sort(sortedData.begin(), sortedData.end());

    //Compute cumulative distribution
    std::vector<float> cdf(sortedData.size(), 0.f);
    cdf[0] = sortedData[0].second;
    for (size_t i = 1; i < sortedData.size(); ++i) {
        cdf[i] = cdf[i - 1] + sortedData[i].second;
    }

    //Determine quantiles
    size_t dataIndex = 0;
    for (const auto& [key, targetCdf] : quantileKeys) {
        while (dataIndex < sortedData.size() && cdf[dataIndex] < targetCdf) {
            ++dataIndex;
        }
        sumStat[key] = sortedData[dataIndex].first;
    }

    return sumStat;
}

Vecteur<float> Admixture::firstMoments(int N)
{
    Vecteur<float> firstMoments(N,0);
    float weightSum = 0.f;
    float weightSumSq = 0.f;

    for (auto j = (this->fraction_H1).begin(); j != (this->fraction_H1).end(); j++)
    {
        weightSum += j->second * j->first;
        weightSumSq += j->second * j->first * j->first;
    }

    firstMoments[0] = weightSum;
    firstMoments[1] = weightSumSq - weightSum*weightSum;
    
    for (int i=2; i<N; i++)
    {   
        for (auto j = (this->fraction_H1).begin(); j != (this->fraction_H1).end(); j++)
        {
            float diff = j->first - weightSum;
            firstMoments[i] += j->second * pow(diff,i+1);
        }

        firstMoments[i] = firstMoments[i] / pow(sqrt(firstMoments[1]), i+1);
    }

    return firstMoments;
}

Admixture& Admixture::RandomMating(Admixture &H_t_tilde_plus_1, Vecteur<float> contribution, Vecteur<Vecteur<float>> e)
{   
    //Number of source population
    int n = contribution.size()-1;
 
    //Contribution-Contribution descendants and Contribution-Hybrid descendants 
    for (int k=0; k<n; k++)
    {
        if (contribution[k]!=0)         //If contribution of population k is null, it isn't involve in admixture at that generation
        {
            for (int l=0; l<n; l++)
            {
                if (contribution[l]!=0) //If contribution of population l is null, it isn't involve in admixture at that generation
                {
                    H_t_tilde_plus_1.fraction[0.5f*(e[k]+e[l])] += contribution[k]*contribution[l];
                }
            } 

            for (const auto& h1 : (*this).fraction)
            {
                H_t_tilde_plus_1.fraction[0.5f*(e[k]+h1.first)] += 2.f*contribution[k]*h1.second*contribution.back();
            }   
        }
    }

    //Hybrid-Hybrid descendants
    for (const auto& parent1 : (*this).fraction)
    {
        for (const auto& parent2 : (*this).fraction)
        {
            if (parent1.second * parent2.second != 0)
            {
                H_t_tilde_plus_1.fraction[0.5f * (parent1.first+parent2.first)] += parent1.second * parent2.second * pow(contribution.back(),2);
            }
        }
    }

    return(H_t_tilde_plus_1);
}

Admixture& Admixture::Selection(Admixture &H_t, Vecteur<Vecteur<float>> Phi)
{
    //Store Phi_t_ro
    map<Vecteur<float>, double> phi_t_ro;

    //Compute alpha_t
    float alpha_t_inv = 0.f;

    //Go through all genotypes and compute phi(t,ro)
    for (const auto& fraction : (*this).fraction)
    {
        int n = (fraction.first).size(); 
        float phi_value = 0.f;

        for (int i=0; i<n; i++)
        {
            for (int j=0; j<n; j++)
            {   
                phi_value += Phi[i][j] * fraction.first[i] * fraction.first[j];
            }
        }

        phi_t_ro[fraction.first] = phi_value;
        alpha_t_inv +=  phi_value * fraction.second;
    }
    //Compute H_t from H_t_tilde and previously calculated alpha_t and Phi_ro_t
    if (alpha_t_inv  == 0.f)
    {
        cout << "error alpha_t_inv is null" << endl;
        exit(1);
    }

    else 
    {
        for (const auto& fraction : (*this).fraction)
        {
            H_t.fraction[fraction.first] = 1.f / alpha_t_inv * phi_t_ro[fraction.first] * fraction.second; 
        }
    }

    return(H_t);
}

Admixture& Admixture::RandomMating_H1(Admixture &H_t_tilde_plus_1, Vecteur<float> contribution)
{
    //Contribution-Contribution descendants and Contribution-Hybrid descendants 
    if (contribution[0]!=0 and contribution[1]!=0)
    {
        H_t_tilde_plus_1.fraction_H1[1.f] = contribution[0]*contribution[0];
        H_t_tilde_plus_1.fraction_H1[0.5f] = 2.f*contribution[0]*contribution[1];
        H_t_tilde_plus_1.fraction_H1[0.f] = contribution[1]*contribution[1];

        for (const auto& h1 : (*this).fraction_H1)
        {
            H_t_tilde_plus_1.fraction_H1[0.5f*(1+h1.first)] += 2.f*contribution[0]*h1.second*contribution.back();
            H_t_tilde_plus_1.fraction_H1[0.5f*(h1.first)] += 2.f*contribution[1]*h1.second*contribution.back();
        }
    }
    
    else if (contribution[0]==0 and contribution[1]!=0)
    {
        H_t_tilde_plus_1.fraction_H1[0.f] = contribution[1]*contribution[1];

        for (const auto& h1 : (*this).fraction_H1)
        {
            H_t_tilde_plus_1.fraction_H1[0.5f*(h1.first)] += 2.f*contribution[1]*h1.second*contribution.back();
        }
    }

    else if (contribution[0]!=0 and contribution[1]==0)
    {
        H_t_tilde_plus_1.fraction_H1[1.f] = contribution[0]*contribution[0];
        
        for (const auto& h1 : (*this).fraction_H1)
        {
            H_t_tilde_plus_1.fraction_H1[0.5f*(1+h1.first)] += 2.f*contribution[0]*h1.second*contribution.back();
        }
    }

    //Hybrid-Hybrid descendants
    for (const auto& parent1 : (*this).fraction_H1)
    {
        for (const auto& parent2 : (*this).fraction_H1)
        {
            if (parent1.second * parent2.second != 0)
            {
                H_t_tilde_plus_1.fraction_H1[0.5f * (parent1.first+parent2.first)] += parent1.second * parent2.second * pow(contribution.back(),2);
            }
        }
    }

    return(H_t_tilde_plus_1);
}

Admixture& Admixture::RandomMating_H1_discrete(Admixture &H_t_tilde_plus_1, Vecteur<float> contribution, Vecteur<float> discretization)
{
    int N = discretization.size();

    //Contribution-Contribution descendants and Contribution-Hybrid descendants 
    if (contribution[0]!=0 and contribution[1]!=0)
    {
        int ind = round((0.5f - discretization[0]) / (discretization[N - 1] - discretization[0]) * (N - 1));
        H_t_tilde_plus_1.fraction_H1[discretization[N-1]] = contribution[0]*contribution[0];
        H_t_tilde_plus_1.fraction_H1[discretization[ind]] = 2.f*contribution[0]*contribution[1];
        H_t_tilde_plus_1.fraction_H1[discretization[0]] = contribution[1]*contribution[1];

        for (const auto& h1 : (*this).fraction_H1)
        {
            int ind1 = round((0.5f*(1+h1.first) - discretization[0]) / (discretization[N - 1] - discretization[0]) * (N - 1));
            int ind2 = round((0.5f*(h1.first) - discretization[0]) / (discretization[N - 1] - discretization[0]) * (N - 1));
            H_t_tilde_plus_1.fraction_H1[discretization[ind1]] += 2.f*contribution[0]*h1.second*contribution.back();
            H_t_tilde_plus_1.fraction_H1[discretization[ind2]] += 2.f*contribution[1]*h1.second*contribution.back();
        }
    }
    
    else if (contribution[0]==0 and contribution[1]!=0)
    {
        H_t_tilde_plus_1.fraction_H1[discretization[0]] = contribution[1]*contribution[1];

        for (const auto& h1 : (*this).fraction_H1)
        {
            int ind = round((0.5f*(h1.first) - discretization[0]) / (discretization[N - 1] - discretization[0]) * (N - 1));
            H_t_tilde_plus_1.fraction_H1[discretization[ind]] += 2.f*contribution[1]*h1.second*contribution.back();
        }
    }

    else if (contribution[0]!=0 and contribution[1]==0)
    {
        H_t_tilde_plus_1.fraction_H1[1.f] = contribution[0]*contribution[0];
        
        for (const auto& h1 : (*this).fraction_H1)
        {
            int ind = round((0.5f*(1+h1.first) - discretization[0]) / (discretization[N - 1] - discretization[0]) * (N - 1));
            H_t_tilde_plus_1.fraction_H1[discretization[ind]] += 2.f*contribution[0]*h1.second*contribution.back();
        }
    }

    //Hybrid-Hybrid descendants
    for (const auto& parent1 : (*this).fraction_H1)
    {
        for (const auto& parent2 : (*this).fraction_H1)
        {
            if (parent1.second * parent2.second != 0)
            {
                int ind = round((0.5f * (parent1.first+parent2.first) - discretization[0]) / (discretization[N - 1] - discretization[0]) * (N - 1));
                H_t_tilde_plus_1.fraction_H1[discretization[ind]] += parent1.second * parent2.second * pow(contribution.back(),2);
            }
        }
    }

    return(H_t_tilde_plus_1);
}

Admixture& Admixture::Selection_H1(Admixture &H_t ,Vecteur<Vecteur<float>> Phi)
{
    //Store Phi_t_ro
    map<float, double> phi_t_ro;

    //Compute alpha_t
    float alpha_t_inv = 0.f;

    //Go through all genotypes and compute phi(t,ro)
    for (const auto& genotype : (*this).fraction_H1)
    {
        float phi_value = 0.f;

        phi_value += Phi[0][0] * genotype.first * genotype.first;
        phi_value += 2.f * Phi[0][1] * genotype.first * (1-genotype.first);
        phi_value += Phi[1][1] * (1-genotype.first) * (1-genotype.first);

        phi_t_ro[genotype.first] = phi_value;
        alpha_t_inv +=  phi_value * genotype.second;
    }

    //Compute H_t from H_t_tilde and previously calculated alpha_t and Phi_ro_t
    if (alpha_t_inv  == 0.f)
    {
        cout << "error alpha_t_inv is null" << endl;
        exit(1);
    }

    else 
    {
        for (const auto& genotype : (*this).fraction_H1)
        {
            H_t.fraction_H1[genotype.first] = 1.f / alpha_t_inv * phi_t_ro[genotype.first] * genotype.second; 
        }
    }

    return(H_t);
}

Admixture Admixture::convert_fraction_H1_to_fraction()
{
    Admixture converted(*this);

    for (const auto& fraction : (*this).fraction_H1)
    {
        converted.fraction[Vecteur<float>({fraction.first})] = fraction.second;
    }

    return(converted);
}

//===========================================================================
//                          External functions
//===========================================================================

ostream &operator<<(ostream &out, const Admixture &g)
{   
    if (g.fraction.size()!=0)
    {
        for (const auto& pair : g.fraction) {
        out << "(" << pair.first << ", " << pair.second << "), ";
        }
    }
    else
    {
        for (const auto& pair : g.fraction_H1) {
            out << "(" << pair.first << ", " << pair.second << "), ";
        }
    }
  return (out);
};

#endif