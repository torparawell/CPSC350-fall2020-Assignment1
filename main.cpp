#include "dna.h"
#include<iostream>
#include<cmath>
//#include<cstdlib>
//#include<math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>

using namespace std;

int main(int argc, char** argv) {
  Calculate *c = new Calculate(); //assigning Calculate to the heap
  fstream newfile;
  bool playAgain = true;
  while(playAgain == true) { //Creates a while loop to execute code.
    //Declaring Variables
    double totalLength = 0.0;
    double tempLength = 0.0;
    double lineNumber = 0.0;
    double letterAcount = 0.0;
    double letterCcount = 0.0;
    double letterGcount = 0.0;
    double letterTcount = 0.0;
    double letterCount = 0.0;
    double mean = 0.0;
    double variance = 0.0;
    double tempVariance = 0.0;
    string currentLine = "";
    float g = 0.00f;
    float d = 0.00f;

    string newDNAstring = "";
    int countBigramAA = 0;
    int countBigramAC = 0;
    int countBigramAG = 0;
    int countBigramAT = 0;
    int countBigramCA = 0;
    int countBigramCC = 0;
    int countBigramCG = 0;
    int countBigramCT = 0;
    int countBigramGA = 0;
    int countBigramGC = 0;
    int countBigramGG = 0;
    int countBigramGT = 0;
    int countBigramTA = 0;
    int countBigramTC = 0;
    int countBigramTG = 0;
    int countBigramTT = 0;
    int totalBigramCount = 0;
    string nucleotides = "";
    string filename = "";
    //Prompting the user for file name.
    cout << "Enter file name: " << endl;
    cin >> filename;
    newfile.open(filename,ios::in);
    if (newfile.is_open()){ //Opening file.
       string tp;
       while(getline(newfile, currentLine)){ //Providing all the contents of the file.
         lineNumber += 1;

         for (int i = 0; i < currentLine.size(); ++i) { //Counting all letters and bigrams.
           if (toupper(currentLine[i]) == 'A') {
             letterAcount++;
           }
           if (toupper(currentLine[i]) == 'C') {
             letterCcount++;
           }
           if (toupper(currentLine[i]) == 'G') {
             letterGcount++;
           }
           if (toupper(currentLine[i]) == 'T') {
             letterTcount++;
           }
         }
         string currentBigram;
         string a1;
         string a2;
         for (int i = 0; i < currentLine.size(); i = i + 2) {
           //Calculating bigram relative frequency.
           if(currentLine.size() - 1 == i) { //Seeing if string length is odd.
             if(currentLine.size() % 2 == 1) {
               a1 = toupper(currentLine[0]);
               a2 = toupper(currentLine.back());
               currentBigram = a1 + a2;
             }
           } else {
             a1 = toupper(currentLine[i]); //Concatenating two chars together to create bigram.
             a2 = toupper(currentLine[i+1]);

             currentBigram = a1 + a2;
           }
           if (currentBigram == "AA") { //If it equals a certain bigram than add to countBigramAA.
             countBigramAA++;
           }
           if (currentBigram == "AC") {
             countBigramAC++;
           }
           if (currentBigram == "AG") {
             countBigramAG++;
           }
           if (currentBigram == "AT") {
             countBigramAT++;
           }
           if (currentBigram == "CA") {
             countBigramCA++;
           }
           if (currentBigram == "CC") {
             countBigramCC++;
           }
           if (currentBigram == "CG") {
             countBigramCG++;
           }
           if (currentBigram == "CT") {
             countBigramCT++;
           }
           if (currentBigram == "GA") {
             countBigramGA++;
           }
           if (currentBigram == "GC") {
             countBigramGC++;
           }
           if (currentBigram == "GG") {
             countBigramGG++;
           }
           if (currentBigram == "GT") {
             countBigramGT++;
           }
           if (currentBigram == "TA") {
             countBigramTA++;
           }
           if (currentBigram == "TC") {
             countBigramTC++;
           }
           if (currentBigram == "TG") {
             countBigramTG++;
           }
           if (currentBigram == "TT") {
             countBigramTT++;
           }
           totalBigramCount++;
           //cout << "Current i: " << i << "  Curent a1: " << currentLine[i] << "  Current a2: " << currentLine[i+1] << endl;
           //cout << "Bigram: " << currentBigram << endl;
         }
         tempLength = c->CalculateLength(currentLine); //Calls the calculate length method.
         totalLength = totalLength + tempLength;
       }
    newfile.close(); //close the file object.
    }

    //Final Information
    //Relativity and frequency
    letterCount = letterAcount + letterCcount + letterGcount + letterTcount;
    double relativeProbabilityLetterA = letterAcount / letterCount; //Calculates relative probabiity of each letter count.
    double relativeProbabilityLetterC = letterCcount / letterCount;
    double relativeProbabilityLetterG = letterGcount / letterCount;
    double relativeProbabilityLetterT = letterTcount / letterCount;

    //Mean of DNA Strings
    mean = totalLength / lineNumber;
    //Variance
    newfile.open("dnaData.txt",ios::in);
    if (newfile.is_open()){
       string tp;
       while(getline(newfile, tp)){
         string currentLine = tp;
         tempLength = c->CalculateLength(currentLine);
         tempVariance += pow((tempLength - mean), 2);
       }
    tempVariance /= lineNumber;
    variance = tempVariance;

    newfile.close(); //close the file object.
    }
    //Standard Deviation
    double standardDeviation = sqrt(variance);

    //Output to out file.
    newfile.open("TorParawell.out",ios::app);
    if(newfile.is_open()) {
      newfile<<"The sum of the length of the DNA strings is: " << totalLength << endl;
      newfile<<"The mean of the length of the DNA strings is: " << mean << endl;
      newfile<<"The variance of the length of the DNA strings is: " << variance << endl;
      newfile<<"The standard deviation of the length of the DNA strings is: " << standardDeviation << endl;
      newfile<<"Here are the relative probability of each nucleotide:" << endl;
      newfile<<"A: " << relativeProbabilityLetterA << endl;
      newfile<<"C: " << relativeProbabilityLetterC << endl;
      newfile<<"G: " << relativeProbabilityLetterG << endl;
      newfile<<"T: " << relativeProbabilityLetterT << endl;
      newfile<<"Here are the relative probabilities of each nucleotide:" << endl;
      newfile<<"AA: " <<(double(countBigramAA) / double(totalBigramCount))<< endl;
      newfile<<"AC: " <<(double(countBigramAC) / double(totalBigramCount))<< endl;
      newfile<<"AG: " <<(double(countBigramAG) / double(totalBigramCount))<< endl;
      newfile<<"AT: " <<(double(countBigramAT) / double(totalBigramCount))<< endl;
      newfile<<"CA: " <<(double(countBigramCA) / double(totalBigramCount))<< endl;
      newfile<<"CC: " <<(double(countBigramCC) / double(totalBigramCount))<< endl;
      newfile<<"CG: " <<(double(countBigramCG) / double(totalBigramCount))<< endl;
      newfile<<"CT: " <<(double(countBigramCT) / double(totalBigramCount))<< endl;
      newfile<<"GA: " <<(double(countBigramGA) / double(totalBigramCount))<< endl;
      newfile<<"GC: " <<(double(countBigramGC) / double(totalBigramCount))<< endl;
      newfile<<"GG: " <<(double(countBigramGG) / double(totalBigramCount))<< endl;
      newfile<<"GT: " <<(double(countBigramGT) / double(totalBigramCount))<< endl;
      newfile<<"TA: " <<(double(countBigramTA) / double(totalBigramCount))<< endl;
      newfile<<"TC: " <<(double(countBigramTC) / double(totalBigramCount))<< endl;
      newfile<<"TG: " <<(double(countBigramTG) / double(totalBigramCount))<< endl;
      newfile<<"TT: " <<(double(countBigramTT) / double(totalBigramCount))<< endl;
      newfile<<"Generated DNA Strings: "<< endl;
    }
    //Creates random function based on computer's time.
    srand(time(NULL));
    //Appends to file 1000 DNA strings.
    for (int l = 0; l < 1001; ++l) {
      newDNAstring = "";
      double a = rand()/double(RAND_MAX);
      double b = rand()/double(RAND_MAX);
      g = sqrt(-2 * log(a)) * cos(2 * atan(1) * 4 * b);
      d = int(standardDeviation * g + mean);
      //Calculates how many letters should be in the string based on relative frequency.
      int lengthA = relativeProbabilityLetterA * d;
      int lengthC = relativeProbabilityLetterC * d;
      int lengthG = relativeProbabilityLetterG * d;
      int lengthT = relativeProbabilityLetterT * d;
      nucleotides = "";
      //Appends the amount to a new string called nucleotides.
      for (int i = 0; i < lengthA; ++i) {
        nucleotides = nucleotides + "A";
      }
      for (int i = 0; i < lengthC; ++i) {
        nucleotides = nucleotides + "C";
      }
      for (int i = 0; i < lengthG; ++i) {
        nucleotides = nucleotides + "G";
      }
      for (int i = 0; i < lengthT; ++i) {
        nucleotides = nucleotides + "T";
      }
      newDNAstring = "";
      //cout << "Nucleotides: " << nucleotides << endl;
      for (int i = 0; i < int(d); ++i) { //Randomizes and scrambles the string using swap method.
        int randNum = rand() % int(d);
        swap(nucleotides[i], nucleotides[randNum]);
      }
      newDNAstring = nucleotides;
      //cout << "New String: " << newDNAstring << endl;
      if (newDNAstring != " ") {
        if (newDNAstring != "") {
          newfile << newDNAstring << endl;
        }
      }
      //cout << "Count: " << l << endl;
    }

    string choice = ""; //Asks the user if they want to play again.
    cout << "Do you want to play again?: (Y/N)" << endl;
    cin >> choice;

    if (choice == "Y") {
      playAgain = true;
      newfile << "\nFile: " << filename << endl; //Shows the user that there is a new file.
    } else {
      playAgain = false;
      newfile.close(); //Closes the file when the user says anything other than "Y"
    }
  }

  delete c;
  return 0;
}
