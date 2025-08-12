#include <random>
#include <vector>
#include <complex>
#include <iostream>



std::random_device rd;
std::uniform_int_distribution<int> genValue(0, 0x7fffffff);

double Rand() {
    return ((double)genValue(rd) / RAND_MAX) * 2 * M_PI;   
}

void convert(std::vector<std::complex<double>> &csig, const std::vector<double> &sig) {
    csig.resize(sig.size());
    for(size_t i = 0 ; i < sig.size(); i++) {
        csig[i] = std::complex<double>(sig[i], 0.0);
    }
}

void genSig(std::vector<double> &sig,
            const double duration, const int rate, const double amp, const double noise_level) {
    double size = rate * duration;
    sig.resize(size);
    
    for(size_t i = 0 ; i < size; i++) {
        const double noise = Rand() * noise_level;
        sig[i] = amp * sin(2 * M_PI * 400 * i / rate) 
               + amp * 0.5 * sin(2 * M_PI * 500 * i / rate)
               + amp * 0.3 * sin(2 * M_PI * 600 * i / rate)
               + amp * 0.1 * sin(2 * M_PI * 700 * i / rate)
               + noise;
    }
}


void DFT(std::vector<std::complex<double>>& csig, const std::vector<double> &sig) {
    csig.resize(sig.size());
    
    for(size_t i = 0 ; i < csig.size(); i++) {
        std::complex<double> sum (0.0, 0.0);
        const double var = 2 * M_PI * i / sig.size();
        for(size_t j = 0 ; j < sig.size(); j++) {
            const double angle = var * j;
            std::complex<double> exp(cos(angle), sin(angle));
            sum += sig[j] * exp;
        }
        csig[i] = sum;
    }
}
std::vector<std::complex<double>> FFT(std::vector<std::complex<double>>& csig) {
    size_t n = csig.size();
    if(n <= 1) return csig; 
    
    std::vector<std::complex<double>> even(n / 2);
    std::vector<std::complex<double>> odd(n / 2);
    
    for(size_t i = 0; i < n / 2; i++) {
        even[i] = csig[i * 2];
        odd[i]  = csig[i * 2 + 1];
    }
    
    even = FFT(even);
    odd  = FFT(odd);
    
    std::vector<std::complex<double>> result(n);
    for(size_t i = 0 ; i < n / 2; i++) {
        double angle = 2 * M_PI * i / n;
        std::complex<double> exp(cos(angle), sin(angle));
        result[i]         = even[i] + exp * odd[i];
        result[i + n / 2] = even[i] - exp * odd[i];
    }

    return result;
}


void printMagnitude(std::vector<std::complex<double>>& csig, const int rate) {
    for(size_t i = 0 ; i < csig.size() / 2 ; i++) {
        double magn = abs(csig[i]) / csig.size();
        if(magn > 0.005)
            std::cout << (double(i) * rate) / csig.size() << "Hz, " << magn << "\n";
    }
}


std::vector<std::complex<double>> ifft(std::vector<std::complex<double>> tsig) {
    size_t n = tsig.size();
    
    if (n == 1) return tsig;
    
    std::vector<std::complex<double>> even(n/2), odd(n/2);
    for (size_t i = 0; i < n/2; i++) {
        even[i] = tsig[2*i];
        odd[i]  = tsig[2*i + 1];
    }
    
    even = ifft(even);
    odd  = ifft(odd);
    
    std::vector<std::complex<double>> result(n);
    for (size_t i = 0; i < n/2; i++) {
        double angle = -2 * M_PI * i / n;
        std::complex<double> exp(cos(angle), sin(angle));
        result[i]         = even[i] + exp * odd[i];
        result[i + n / 2] = even[i] - exp * odd[i];
    }
    
    return result;
}
std::vector<std::complex<double>> IFFT(const std::vector<std::complex<double>>& tsig) {
    std::vector<std::complex<double>> result = ifft(tsig);

    for (auto& i : result) {
        i /= tsig.size();
    }
    return result;
}


int main(int argc, char **argv) {
    const int rate = 1024 * 64;

    std::vector<double> sig;
    std::vector<std::complex<double>> tsig;
    
    genSig(sig, 8, rate, .5, 0.01);
    
    
    for(int i = 0; i < 1000; i++) std::cout << sig[i] << " ";
    std::cout << "\n\n";

    convert(tsig, sig);
    tsig = FFT(tsig);
    
    printMagnitude(tsig, rate);
    
    
    std::vector<std::complex<double>>itsig = IFFT(tsig);
    
    for(int i = 0; i < 1000; i++) std::cout << itsig[i].real() << " ";
    std::cout << "\n\n";
    
    
    std::cout << sig.size() << "\n";
    std::cout << itsig.size() << "\n\n";
    
    
    for(size_t i = 0; i < sig.size(); i++)
        std::cout << sig[i] - itsig[i].real() << " ";
    std::cout << "\n\n";
    
    
    return 0;
}
