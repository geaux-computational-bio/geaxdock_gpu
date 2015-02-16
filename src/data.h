/*
==============================================================================================
     __________________ ____ ___              .___             __      _________   _____   
    /  _____/\______   \    |   \           __| _/____   ____ |  | __ /   _____/  /     \  
   /   \  ___ |     ___/    |   /  ______  / __ |/  _ \_/ ___\|  |/ / \_____  \  /  \ /  \ 
   \    \_\  \|    |   |    |  /  /_____/ / /_/ (  <_> )  \___|    <  /        \/    Y    \
    \______  /|____|   |______/           \____ |\____/ \___  >__|_ \/_______  /\____|__  /
           \/                                  \/           \/     \/        \/         \/ 

      GPU-accelerated hybrid-resolution ligand docking using Replica Exchange Monte Carlo

==============================================================================================
*/


#ifndef __DATA_H_
#define __DATA_H_

#include<string>
#include<iostream>
#include<cstdlib>

std::string three2oneS (std::string);

char three2oneC (std::string);

std::string one2three (std::string);

int getResCode (std::string);

int getResCodeOne (std::string);

std::string getResName (int);

int getPntCode (std::string);

std::string getPntName (int);

int getLigCode (std::string);

std::string getLigName (int);

#endif
