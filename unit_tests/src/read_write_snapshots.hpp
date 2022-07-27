#include <ostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>
/**
 * @brief Object to handle snapshot data reading and writing. 
 * A snapshot is static data to compare to during test cases.
 * Data is always loaded into a std::vector class for ease of use.
 * 
 * @tparam T datatype of snapshot data, could be for example int, float, string
 */
template<typename T>
class Snapshot {
private:
    std::vector<T> data;
    std::string filename;
public:
    Snapshot(std::vector<T> data, std::string filename):
        data(data), filename(filename) {}

    Snapshot(std::string filename):
        filename(filename) {}

    void write_to_file() {
        std::ofstream outfile(filename);
        std::ostream_iterator<T> out_iterator(outfile,"\n");
        std::cout << filename << '\n';
        copy(data.begin(),data.end(),out_iterator);
        //std::cout << "here \n";
    };

    std::vector<T> read_from_file() {
        std::ifstream fin(filename);
        std::istream_iterator<int> fin_it(fin);
        std::istream_iterator<int> eos;
        copy(fin_it, eos, back_inserter(data));
        return data;
    };

    void clear_data() {
        data.clear();
    }
};

