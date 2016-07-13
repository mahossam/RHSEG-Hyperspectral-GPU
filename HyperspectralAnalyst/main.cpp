#include <QtGui>
#include <QApplication>


#include "myglobals.h"
#include "CPUAlgorithm.h"

#include "Libraries\HSWO.h"

#if USE_AMP
#include <amp_math.h>
#include <amp.h> 
//#include <amp_stl_algorithms.h>
#endif 

#include  <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	CPUAlgorithm form;	
    form.show();    
	printf("hi\n");
	
	#if USE_AMP
	auto acc = Concurrency::accelerator();
	wcout<< "get_description "<< acc.get_description() <<endl;
	wcout<< "get_device_path "<<acc.get_device_path() <<endl;
	cout<< "get_version "<<acc.get_version() <<endl;	
	cout<< "get_is_emulated "<<acc.get_is_emulated() <<endl;	
	cout<< "get_dedicated_memory "<<acc.get_dedicated_memory() <<endl;	
	cout<< "get_supports_double_precision "<<acc.get_supports_double_precision() <<endl;
	cout<< "get_supports_limited_double_precision "<<acc.get_supports_limited_double_precision() <<endl;
	#endif 

	if(argc > 1)
	{
		HyperSpectralToolbox::HSWO::silentMode = true;
	}

	//Exiting CUDA utility library initiated inside cuda code by cutilDeviceInit ...
	//cutilExit(NULL, NULL);

	return a.exec();
}
