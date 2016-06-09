#include "FluidSimulator.hh"
#include "Types.hh"


int main( int argc, char** argv )
{
    CHECK_MSG(argc==2, "Filename missing");
    std::cout<<"\nReading file " << argv[1] <<std::endl;
	

    FileReader input_data;

    input_data.readFile(argv[1]);

    FluidSimulator simulate(input_data);

    real steps = input_data.getIntParameter("timesteps");
    simulate.simulateTimeStepCount(steps);

}
