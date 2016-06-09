
#include "FileReader.hh"



void FileReader::registerIntParameter(const std::string &key, int init)
{
   std::stringstream ss;
   ss << init;
   parameter[key] = ss.str();
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
   std::stringstream ss;
   ss << init;
   parameter[key] = ss.str();
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
   parameter[key] = init;
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
   CHECK_MSG(parameter.find(key) != parameter.end(), "Parameter '" + key + "\"' is not in registery");
   parameter[key] = in;
}

void FileReader::setParameter(const std::string &key, real in)
{
   std::stringstream ss;
   ss << in;
   CHECK_MSG(parameter.find(key) != parameter.end(), "Parameter '" + key + "\"' is not in registery.");
   parameter[key] = ss.str();
}

void FileReader::setParameter(const std::string &key, int in)
{
   std::stringstream ss;
   ss << in;
   CHECK_MSG(parameter.find(key) != parameter.end(), "Parameter '" + key + "\"' is not in registery.");
   parameter[key] = ss.str();
}


bool FileReader::readFile(const std::string &name)
{
   std::ifstream Data_file(name);
   CHECK_MSG(Data_file.is_open(),"Could not open file '"+ name + "' which has to be in the current directory");
   std::string row, col_0, col_1;
   while (std::getline(Data_file, row)) 
   {
       std::istringstream ss(row, std::istringstream::in);
       ss >> col_0 >> col_1;

       if ( !col_0.empty() && col_0.at(0) != '#' && !col_1.empty() && col_1.at(0) != '#')
       {
          if(parameter.find(col_0) == parameter.end())
             registerStringParameter(col_0,col_1);
          else
          {
             try {stod(parameter[col_0]);}
             catch(...){
                  setParameter(col_0,col_1);
                  continue;}
             try {stod(col_1);
              }
             catch(...)
              {
               CHECK_MSG(0, "Type mismatch for parameter '" + col_0 +"' "); }
             setParameter(col_0,real(stod(col_1))); }
       }
   }   
   PROGRESS("File reading complted .")
   return false;
}



void FileReader::printParameters() const
{
   for(auto it = parameter.begin(); it != parameter.end(); ++it)
    std::cout << std::left<<"\t\t"<< it->first << "      " << it->second << "\n";
}

bool FileReader::checkparameter(const std::string & key ) const
{
   if (parameter.find(key) != parameter.end())
       return true;
   else
       return false;
}
