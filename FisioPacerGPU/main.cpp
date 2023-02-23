/* 
 * File:   main.cpp
 * Author: ricardo
 *
 * Created on October 24, 2014, 7:56 PM
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include <omp.h>

#include "kernel.h"
//#include "Stopwatch.h"



/**
 * 
 * @param argc
 * @param argv
 * @return 
 */
int main(int argc, char** argv) {
//    printf("FIsioPacer2.0\n");getchar();
    using namespace std;
    string paramFile;
    if(argc != 2)
    {
        setDefaultFolders();
        cout << "Correct usage:" << endl;
        cout << "\tfisiopacerrun paramFile" << endl;
        exit(0);
    } else{
        paramFile       = argv[1];
    }
//./dist/Release/GNU-Linux-x86/fisiopacerrun /home/ricardo/tese/cubo6/malha.param
    //Stopwatch watch;
   // watch.start();
    startCA(paramFile, true);
    //watch.stop();
    //cout<<"Time "<< watch.timeMS() <<" "<< (watch.timeMS()/1000)/60 
    //<<"'"<<(watch.timeMS()/1000)%60<<endl;
    return 0;
}

