// This function calculates the highest number of pairs of nitrougenous bases in a DNA helix (AT, GC) of each str in the loaded database.
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "ourvector.h"
using namespace std;

// DNAStore - Struct to store each person with their corresponding STR counts
struct DNAStore {
    string name;
    ourvector<int> DNAvalues;
};

//////// MILESTONE 1
// LineOne_input()
// This function adds STR combinations from the first line of the text file into a container
void LineOne_input(string lineOne, ourvector<ourvector<char>>& gene, ourvector<DNAStore>& personDNA, bool& load) {
    if (load == true) {  // If previously loaded, clear the vectors
        gene.clear();
        personDNA.clear();
    }
    stringstream check(lineOne);
    string name; 
    string num;
    getline(check, name, ',');  // 'name' is avoided

    while (getline(check, num, ','))
    {
      ourvector<char> value;
      for (char i: num) {
        value.push_back(i);
      }
      gene.push_back(value);
    }
}

// load_db()
// This function loads the database into appropriate containers
void load_db(ifstream& inFile, string filename, ourvector<ourvector<char>>& gene,ourvector<DNAStore>& personDNA, bool& load){
    cout<<"Loading database..." << endl;
    inFile.open(filename);
    if(!inFile.is_open()){
        cout<<"Error: unable to open '"<< filename <<"'"<<endl;
        return;
    }
    string lineOne, allLines;  // first line is different to the rest of the lines
    DNAStore profile;
    inFile >> lineOne;
    LineOne_input(lineOne, gene, personDNA, load);  // storing first line into gene
    while (getline(inFile, allLines)) {  // looping through rest of the lines
        stringstream check(allLines);
        string name, stringNum;
        getline(check, name, ',');
        profile.name = name;  // storing name of each person
        profile.DNAvalues.clear();
        while (getline(check, stringNum, ','))
        {
          if (stringNum != ""){  // incase empty
            int num = stoi(stringNum);
            profile.DNAvalues.push_back(num);  // storing amount of dna for each person
          }
        }
        personDNA.push_back(profile);  // adding to overall vector
    }
    inFile.close();
}

//////// MILESTONE 3
// load_dna()
// This function loads the dna file into appropriate containers
void load_dna(ourvector<char>& DNA_store, string filename, bool& dna_load) {
  DNA_store.clear();  // making sure it's clear incase it's called before
  cout << "Loading DNA..." << endl;
  string line;
  ifstream inFile;
  inFile.open(filename);
  if (!inFile.is_open()){
    cout << "Error: unable to open '" << filename << "'" << endl;
    return;
  }
  while(getline(inFile, line)) {
    for (int i = 0; i < line.length(); i++) {
      DNA_store.push_back(line[i]);
    }
  }
  dna_load = true;
}

//////// MILESTONE 4
// compare()
// Compares if sequence matches dna chain in the file
int compare(ourvector<char> &gene, ourvector<char> &DNA_store, int index) {
  int i = 0; 
  while (i < gene.size()) {
    if (index >= DNA_store.size()) {
      return 0;
    }
    if (gene[i] != DNA_store[index]) {
      return 0;
    }
    index++;
    i++;
  }
  return 1;
}

// process()
// Identifying the person a given sequence belongs to
void process(ourvector<ourvector<char>>& gene,ourvector<DNAStore>& personDNA, ourvector<char>& DNA_store, ourvector<int>& gene_code){
  if(personDNA.size() < 1) {
    cout << "No database loaded." << endl;
    return;
  } if (DNA_store.size() < 1) {
    cout << "No DNA loaded." << endl;
    return;
  }
  cout << "Processing DNA..." << endl;
  gene_code.clear();  // clearing container incase used before
  for (int i = 0; i < gene.size(); i++) {
    int num = 0, max = 0;
    for (int j = 0; j < DNA_store.size(); j++) {
      if (DNA_store[j] == gene[i][0]) {
        if (compare(gene[i], DNA_store, j) == 1) {  // finding match
          num++;
          if (j == (DNA_store.size() - gene[i].size())) {
            if (max < num) {  // storing the longest sequence
              max = num;
            }
          }
          int gene_size = gene[i].size() - 1;
          j = j + gene_size;
        }
      } else {
        if (max < num) {
          max = num;
        }
        num = 0;
      }
    }
    gene_code.push_back(max);
  }
}

//////// MILESTONE 2
// display_process()
// Prints out the contents after processing vector
void display_process(ourvector<ourvector<char>>& gene, ourvector<int>& gene_code) {
  for (int i = 0; i < gene.size(); i++) {
      for (int j = 0; j < gene[i].size(); j++) {
            cout << gene[i][j];
      }
      cout << ": " << gene_code[i] << endl; 
    }
}

// printDatabase()
// Prints out the contents of the struct database
void printDatabase(ourvector<DNAStore>& personDNA) {
  for (int i = 0; i < personDNA.size(); i++) {
    cout << personDNA[i].name << " ";
    for (int j = 0; j < personDNA[i].DNAvalues.size(); j++) {
          cout << personDNA[i].DNAvalues[j] << " ";
      }
    cout << endl;
  }
}

// printDNA()
// Prints out the contents of the long sequence of dna
void printDNA(ourvector<char> & DNA_store){
  cout << "DNA loaded:" << endl;
  for (int i = 0; i < DNA_store.size(); i++) {
    cout << DNA_store[i];
  }
}

// display()
// Displays to the screen the database, the dna, and STR count
void display(ourvector<ourvector<char>>& gene,ourvector<DNAStore>& personDNA, bool& dna_load, ourvector<char>& DNA_store, bool processed, ourvector<int>& gene_code) {
  if (personDNA.size() < 1) {  // checking if database has loaded
      cout << "No database loaded." << endl;
    if (dna_load == true) {  // checking if dna has been loaded
      printDNA(DNA_store);
      cout << endl << "No DNA has been processed." << endl;
    } else {
      cout << endl << "No DNA loaded." << endl << "No DNA has been processed." << endl;  // dna is only processed if database and dna have been loaded
    }
  } else {
    cout << "Database loaded:" << endl;
    printDatabase(personDNA);
    if (dna_load == false) {
      cout << "No DNA loaded." << endl;
      cout << "No DNA has been processed." << endl;
    }
    else{
      printDNA(DNA_store);
      cout << endl << endl;
      if (processed == false) {  // checking if dna has been processed
        cout << "No DNA has been processed." << endl;
      } else {
        cout << "DNA processed, STR counts:" << endl;
        display_process(gene, gene_code);
      }
    }
  }
}

//////// MILESTONE 5
// search()
// Used the STR counts found previously to search the database and find a match
void search(ourvector<ourvector<char>> &gene,ourvector<DNAStore> &personDNA, ourvector<char> &DNA_store, ourvector<int> &gene_code) {
  if (personDNA.size() < 1) {  // checks if database is loaded
      cout << "No database loaded." << endl;
  } else {
      if (DNA_store.size() < 1) {  // checks if dna is loaded
          cout << "No DNA loaded." << endl;
      } else {
          if (gene_code.size() < 1) {  // checks if dna has been processed
              cout << "No DNA processed." << endl;
          } else {
            cout << "Searching database..." << endl;
            int compare = 0;
            for(int i = 0; i<personDNA.size(); i++) {
              for(int j = 0; j<personDNA[i].DNAvalues.size(); j++) {
                if(gene_code[j] != personDNA[i].DNAvalues[j]) {  // checks if dna count is not the same processed dna count
                  compare = 0;
                  break;
                }
              }
              if(compare == 1) {
                cout<<"Found in database!  DNA matches: "<< personDNA[i].name << endl;
                return;
              }
              compare = 1;
            }
            cout << "Not found in database." << endl;
          }
        }
    }
}

//////// MILESTONE 6
// creative()
// Prior functions required: load_db
// This function primarly looks at pairs of nitrougenous bases in a DNA helix, Adenine (A) only pairs with Thymine (T) and Cytosine (C) only pairs with Guanine (G), and it will count each pair in each STR of the loaded file.
void creative(ourvector<ourvector<char>>& gene) {
  cout << endl << "***Creative function***" << endl << endl;
  for (int i = 0; i < gene.size(); i++) {   // looping through str container
    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    string s = "";
    for (int j = 0; j < gene[i].size(); j++) {
      if(gene[i][j] == 'A'){  // checking each of the base
        count1++;
      } else if (gene[i][j] == 'T'){
        count2++;
      } else if (gene[i][j] == 'C'){
        count3++;
      } else if (gene[i][j] == 'G'){
        count4++;
      }
      s += gene[i][j];  // storing the str characters in a string
    }
    int low1, low2;
    if (count1 < count2) {  // finding the most pairs in the STR of AT
      low1 = count1;
    } else {
      low1 = count2;
    }
    if (count3 < count4) {  // finding the most pairs in the STR of GC
      low2 = count3;
    } else {
      low2 = count4;
    }
    cout << "For " << s << ": AT pairs = " << low1 << ", GC pairs = " << low2 << endl; 
  }
}

// commands()
// A menu system that the user interacts with to manipulate the container accordingly
void commands(ifstream& inFile, string filename, string input, ourvector<ourvector<char>>& gene,ourvector<DNAStore>& personDNA, ourvector<char>& DNA_store, bool& load, bool& dna_load, bool& processed, ourvector<int> gene_code){
    while(input != "#"){

        if(input == "load_db"){
          cin >> filename;
            load_db(inFile, filename, gene, personDNA, load);
            load = true;
        } else if (input == "display"){
            display(gene, personDNA, dna_load, DNA_store, processed, gene_code);
        } else if (input == "load_dna") {
          cin >> filename;
          load_dna(DNA_store, filename, dna_load);
        } else if (input == "process") {
          process(gene, personDNA, DNA_store, gene_code);
          processed = true;
        } else if (input == "search") {
          search(gene, personDNA, DNA_store, gene_code);
        } else if (input == "creative") {
          creative(gene);
        }
        else{
          cout << "incorrect" << endl;
        }
        cout<<"Enter command or # to exit: ";
        cin >> input;
    }
}


int main() {
  cout << "Welcome to the DNA Profiling Application." << endl;
    
    cout<<"Enter command or # to exit: ";
    string input,filename;
    cin >> input;

    ifstream inFile;
    ourvector<ourvector<char>> gene;  // storing all first line STRs
  
    ourvector<DNAStore> personDNA;  // storing name of people and corresponding str counts
    ourvector<char> DNA_store;  // storing long chain of dna
    ourvector<int> gene_code;  // storing processed dna
    bool load = false, dna_load = false, processed = false;
    commands(inFile, filename, input, gene, personDNA, DNA_store, load, dna_load, processed, gene_code);
    return 0;
}