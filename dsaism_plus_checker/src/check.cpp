#include "common.h"
#include "checker.h"

#include <iostream>
#include <chrono>
#include <fstream>

int main(int argc, char** argv){

	// mount stxxl disk	
	stxxl::config *cfg = stxxl::config::get_instance();

	stxxl::disk_config disk("stxxl.tmp", 1024 * 1024 * 1024, "syscall autogrow unlink");

	disk.direct = stxxl::disk_config::DIRECT_ON;

	cfg->add_disk(disk);
	
	// check if input params are legal
	if (argc != 3) {

		std::cerr << "two param required: input_path and output_path.\n";

		exit(-1);
	}	

	// retrieve file name for input string
	std::string s_fname(argv[1]); 
	
	// retrieve file name for output SA
	std::string sa_fname(argv[2]);
	
	// compute the SA for the given string
	Checker<uint8, uint40> checker(s_fname, sa_fname);	

	std::cerr << (checker.run() ? "right" : "wrong")<< std::endl;
}
