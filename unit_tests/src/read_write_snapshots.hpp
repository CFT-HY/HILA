#include <ostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>

template<typename T>
class Snapshots {
private:
    std::vector<T> data;
    std::string filename;
public:
    Snapshots(std::vector<T> data, std::string filename):
        data(data), filename(filename) {}

    void write_to_file() {
        std::ofstream outfile(filename);
        std::ostream_iterator<T> out_iterator(outfile,"\n");
        std::cout << filename << '\n';
        copy(data.begin(),data.end(),out_iterator);
        //std::cout << "here \n";
    };

    std::vector<int> read_from_file() {
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

