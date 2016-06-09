#ifndef FILEREADER_HH
#define FILEREADER_HH

#include<iostream>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>

#include "Types.hh"
#include "Debug.hh"


//*******************************************************************************************************************
/*! Class for reading configuration from a file
*
* Configuration File Syntax:
*   - everything between a '#' character and the beginning of the next line is ignored
*   - lines can be empty
*   - line contain parameter key followed by white-spaces followed by their value
*
*  All possible keys (and the datatype of the value) have to be registered first:
*   - for example usage have a look at the FileReaderTest
*
*
*
*  This Skeleton is just a suggestion. If you are familiar with template programming
*  a possible alternative would be a version that does not need registering of all keys
*  and has a getter function like the following:
*      template<typename ExpectedType>
*      ExpectedType getValue( const std::string & key );
*/
//*******************************************************************************************************************
class FileReader
{
private:

	std::map<std::string,std::string> parameter;

public:

	//register a new parameter with name key and initial int value
	void registerIntParameter( const std::string & key, int init = 0 );

	//register a new parameter with name key and initial double value
	void registerRealParameter( const std::string & key, real init = 0 );

	//register a new parameter with name key and initial string value
	void registerStringParameter( const std::string & key, const std::string & init = "" );

	//set a value for the key string with value in
	void setParameter( const std::string & key, const std::string & in );

	//set a value for the key string with value in
	void setParameter( const std::string & key, real in );

	//set a value for the key string with value in
	void setParameter( const std::string & key, int in );

	// get the int value of key 
	inline int getIntParameter( const std::string & key ) const;

	// get the double value of key 
	inline real getRealParameter( const std::string & key ) const;

	// get the string value of key 
	inline std::string getStringParameter( const std::string & key ) const;

	//try to read all registered parameters from file name
	bool readFile( const std::string & name );

	//print out all parameters to std:out
	void printParameters() const;

    //find if parameter is present
    bool checkparameter( const std::string & key ) const;

private:

};




inline int FileReader::getIntParameter(const std::string &key) const
{
	// take  pointer to parameter
	auto iterator= parameter.find(key); 
	
	// initializing a dummy value
    int value =0;

	//checks if itertator is at end of file , then parameter is not present in file
	CHECK_MSG(iterator != parameter.end(),"Sorry, parameter '"+key+"' does not exist");

	//conversion of string to int
	try {value = std::stoi(iterator->second);}
	catch (...)
       {CHECK_MSG(0, "Parameter passed in " + key +" is not of type int.");}
	return value;
}


inline real FileReader::getRealParameter(const std::string &key) const
{
	auto iterator = parameter.find(key);
    real value=0;
	CHECK_MSG(iterator != parameter.end(), "Sorry parameter '" + key +"' does not exist!");

	//conversion of string to float...Real number are either float or double
	try {value = std::stof(iterator->second);}
	catch (...)
       {CHECK_MSG(0, "Parameter passed in " + key +" is not of type real.");}
	return value;
}

inline std::string FileReader::getStringParameter(const std::string &key) const
{
	auto iterator = parameter.find(key);
	std::string value;
	CHECK_MSG(iterator != parameter.end(), "Error: Parameter '" + key +"' does not exist!");
	value = iterator->second;
	return value;
}





#endif //FILEREADER_HH

