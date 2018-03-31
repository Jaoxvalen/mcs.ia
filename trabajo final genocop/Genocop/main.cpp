#include "Genocop.h"

#include <iostream>
#include <string>

using namespace std;

Genocop* genocop;
int main(int argc, char** argv)
{
    string input;

    if(argc < 2) {
        input = "t02_2";
    }
    else
    {
        input = string(argv[1]);
    }

    genocop = new Genocop("../" + input);
    return 0;
}