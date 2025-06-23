#include <iostream>
#include "assembly.hpp"

int main(int argc, char **argv)
{
    if (argc < 2) {
	std::cout << "please provide a configuration file" << std::endl;
    } else {
	Assembly *assembly = new Assembly(argv[1]);
	assembly->run();
    }
}
