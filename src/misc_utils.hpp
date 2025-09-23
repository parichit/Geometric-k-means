#include <iostream>
#include <vector>
#include <map>
#pragma once

using namespace std;

class print_utils{
        template <typename T> 
        void print_3d_vector(vector<vector<vector<T> > > &data, int num_records,
        string dataname);

        template <typename Tfloat, typename Tint> 
        void print_neighbors(vector<vector<pair<Tfloat, Tint> > > &data, int num_records,
        string dataname);

        template <typename T>
        void print_2d_vector(vector<vector<T> > &, int, string);
        
        template <typename T> 
        void print_vector(vector<T> &, int, string);

        template <typename T1>
        void print_map(map<T1, vector<T1> > &assign_dict, int num_records, 
        string dataname);
};


template <typename Tfloat, typename Tint>  
void print_neighbors(vector<vector<pair<Tfloat, Tint> > > &data, int num_records,
string dataname){

cout << "Printing: " << dataname << " \n" ;

int limit = 0;

cout << "\n";

for (int i=0; i<data.size(); i++){
    
    if (limit < num_records){
        cout << "Record- " << i << "\t" ;
        
        for (int j=0; j<data[i].size(); j++)
            cout << data[i][j].first << ":" << data[i][j].second << "\t";
        }
    
    else{
        break;
    }
        cout << "\n";
        limit++;
    } 

    cout << "\n";
}



template <typename T> void print_2d_vector(vector<vector<T> > &data, int num_records,
string dataname){

cout << "Printing: " << dataname << " \n" ;

int limit = 0;

cout << "\n";

for (int i =0; i< data.size(); i++){
    if (limit < num_records){
        cout << "Record- " << i << "\t" ;
        for (int j=0; j<data[i].size(); j++)
            cout << data[i][j] << "\t";
        }
    else{
        break;
    }
        cout << "\n";
        limit++;
    } 

    cout << "\n";
}


template <typename T> void print_3d_vector(vector<vector<vector<T> > > &data, int num_records,
string dataname){

cout << "Printing: " << dataname << " \n" ;

int limit = 0;

cout << "\n";

for (int i=0; i<data.size(); i++){
    
    cout << "Printing for record: " << i << endl;
    cout << "------------\n";

    if (limit < num_records){
        for (int j=0; j<data[i].size(); j++){
            cout << i <<"-" << j << "\t" ; 
            for (int k = 0; k<data[i][j].size(); k++){
                cout << data[i][j][k] << "\t";
            }
            cout << "\n"; 
        }
    }
    else{
        break;
    }
        cout << "------------\n";
        limit++;
    }
    cout << "\n"; 
}


template <typename T> void print_vector(vector<T> &data, int num_records, string dataname){

cout << "Printing: " << dataname << " \n" ;
int limit = 0;

cout << "\n";

for (int i =0; i< data.size(); i++){
    if (limit < num_records){
        cout << data[i] << "\n";
        limit++;
    }
    else{
        break;
    }
    }
cout << "\n";
}


template <typename T1>
void print_map(map<T1, vector<T1> > &assign_dict, int num_records, string dataname){

cout << "Printing: " << dataname << " \n" ;
int limit = 0;

for(map<int, vector<int> >::iterator ii=assign_dict.begin(); ii!=assign_dict.end(); ++ii){
       cout << (*ii).first << ": ";
       vector <T1> inVect = (*ii).second;
       for (unsigned j=0; j<inVect.size(); j++){
           cout << inVect[j] << " ";
       }
       cout << endl;
   }
}

  //cout << "Printing first five records" << "\n";
    // int limit = 5;
    // int i=0;
    // for(auto row: dataset){
    //     if (i < limit)
    //         for (auto col: row)
    //             cout << col << "\t" ;
    //     else
    //         break;
    //     cout << labels[i] << "\n";
    //     i++;
    // }