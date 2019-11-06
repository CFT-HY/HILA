
#include<string>
#include<map>
#include<list>

/////////////////
/// Class for parsing runtime parameter files 
/// - TODO: robust get function for values. Perhaps somehow make them avaibale?
/// - TODO: implement construct essentials
/// - TODO: implement check --//--     
/////////////////

class inputs {
    public:

        //read runtime parameters from file
        inputs(std::string & fname);
        ~inputs(){};
    
        //return a reference to parameter stored in this object  

        //add an essential variable - triggers complaint if not in parameter file
        void add_essential(std::string & );
        template< typename T >
        //same as above, except defines a default value to be used in case it is not found
        void add_essential(std::string &, T default_value);

    private:

        std::map<std::string, bool> essentials;
        std::map<std::string, double> values;
        std::map<std::string, std::string> names;
        //check that essential variables were found
        void define_essentials(){};
        //handle one line of input
        void handle(std::string &);
        //automatically define necessary runtime variables
        void construct_essentials(){};
};
