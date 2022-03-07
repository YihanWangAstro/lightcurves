#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <thread>
constexpr double pi = 3.14159265358979323846;
constexpr double c_cgs = 3e10;
constexpr double compton_wave_len = 2.43e-10;

inline static double Uniform(double low, double high) {
    static thread_local std::mt19937 generator{std::random_device{}()};
    std::uniform_real_distribution<double> dist{low, high};
    return dist(generator);
}

struct Disk {
   public:
    Disk(double H, double h, double n0) : H_max{H}, scale_height{h}, central_den{n0} {}
    bool in_disk(double x, double y, double z) const {
        if ((fabs(z) >= H_max)) {
            return false;
        } else {
            return true;
        }
    }
    double get_num_den(double x, double y, double z) const {
        return central_den * exp(-z * z / (2 * scale_height * scale_height));
    }

    double H_max{1};
    double scale_height{1};
    double central_den{1};
};

struct Photon {
   public:
    double x{0}, y{0}, z{0};
    double vx{0}, vy{0}, vz{0};
    double t0{0};
    double nu{0};
    double weight{0};
};

double get_mean_free_path(Disk const& disk, Photon& p, double cross_section) {
    double n = disk.get_num_den(p.x, p.y, p.z);
    return 1 / (n * cross_section);
}

double calc_travel_path(double mean_free_path) { return -mean_free_path * log(Uniform(0, 1)); }

double scattering(Disk const& disk, Photon& p, double cross_section) {
    double mfp = get_mean_free_path(disk, p, cross_section);
    double dr = calc_travel_path(mfp);

    double phi = Uniform(0, 2 * pi);

    double cos_theta = Uniform(-1.0, 1.0);
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    if (dr > 2 * disk.H_max / fabs(cos_theta)) {
        dr = 2 * disk.H_max / fabs(cos_theta);
    }

    double dz = dr * cos_theta;
    double dx = dr * sin_theta * cos_phi;
    double dy = dr * sin_theta * sin_phi;

    double v_p = sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);

    p.x += dx, p.y += dy, p.z += dz;

    double new_vx = v_p * sin_theta * cos_phi;
    double new_vy = v_p * sin_theta * sin_phi;
    double new_vz = v_p * cos_theta;

    /* double cos_scatter_angle = (p.vx * new_vx + p.vy * new_vy + p.vz * new_vz) / (v_p * v_p);

     double wave_len_shift = compton_wave_len * (1 - cos_scatter_angle);

     p.nu = 1.0 / (wave_len_shift / c_cgs + 1 / p.nu);*/

    p.vx = new_vx, p.vy = new_vy, p.vz = new_vz;

    p.t0 += dr / c_cgs;

    return mfp;
}

size_t load_photons(std::string const& fname, std::vector<Photon>& photons) {
    std::ifstream file(fname, std::ifstream::in);
    size_t num_line = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
    file.seekg(0);

    std::cout << "proceeding " << num_line << " photons\n";
    photons.reserve(num_line);
    char end_of_line;
    for (size_t i = 0; !file.eof(); i++) {
        Photon p;
        file >> p.t0 >> p.x >> p.y >> p.z >> p.vx >> p.vy >> p.vz /*>> p.nu*/ >> p.weight >> end_of_line;
        p.vx *= c_cgs;
        p.vy *= c_cgs;
        p.vz *= c_cgs;

        photons.emplace_back(p);
    }
    return num_line;
}

void trace_back_to_surface(Disk const& disk, Photon& p) {
    double dt = (fabs(p.z) - disk.H_max) / fabs(p.vz);
    p.x -= p.vx * dt;
    p.y -= p.vy * dt;
    p.z -= p.vz * dt;
    p.t0 -= dt;
}

bool in_vision(double theta, double view_angle, double open_angle) {
    return (theta < view_angle + 0.5 * open_angle) & (theta >= view_angle - 0.5 * open_angle);
}

void evolves(size_t job_id, std::vector<Photon>& photons, Disk const& disk, double cross_section, size_t num_collected,
             double t_max, double view_angle, double open_angle) {
    std::ofstream file("output" + std::to_string(view_angle) + " " + std::to_string(job_id) + ".txt",
                       std::ofstream::out);
    file << std::setprecision(16);
    for (size_t i = 0; i < num_collected;) {
        int j = int(Uniform(0, photons.size()));
        auto& p = photons[j];
        for (;;) {
            double mpf = scattering(disk, p, cross_section);
            if (!disk.in_disk(p.x, p.y, p.z)) {
                trace_back_to_surface(disk, p);

                double theta = acos(p.vz / c_cgs);

                if (in_vision(theta, view_angle * pi / 180.0, open_angle * pi / 180.0)) {
                    file << i << ' ' << p.t0 << ' ' << p.x << ' ' << p.y << ' ' << p.z << ' ' << p.vx / c_cgs << ' '
                         << p.vy / c_cgs << ' ' << p.vz / c_cgs << ' ' << p.nu << ' ' << mpf << ' ' << p.weight << '\n';
                    i++;
                }
                break;
            } else if (p.t0 > t_max) {
                break;
            }
        }
    }
}

void full_evolves(size_t job_id, std::vector<Photon>& photons, Disk const& disk, double cross_section, size_t start_id,
                  size_t end_id, double t_max) {
    std::ofstream file("11full-" + std::to_string(job_id) + ".txt", std::ofstream::out);
    file << std::setprecision(16);
    for (size_t i = start_id; i < end_id; i++) {
        // int j = int(Uniform(0, photons.size()));
        int j = i;
        auto& p = photons[j];
        for (;;) {
            double mpf = scattering(disk, p, cross_section);
            if (!disk.in_disk(p.x, p.y, p.z)) {
                trace_back_to_surface(disk, p);

                file << i << ' ' << p.t0 << ' ' << p.x << ' ' << p.y << ' ' << p.z << ' ' << p.vx / c_cgs << ' '
                     << p.vy / c_cgs << ' ' << p.vz / c_cgs << ' ' << p.nu << ' ' << mpf << ' ' << p.weight << '\n';

                break;
            } else if (p.t0 > t_max) {
                file << i << ' ' << p.t0 << ' ' << p.x << ' ' << p.y << ' ' << p.z << ' ' << p.vx / c_cgs << ' '
                     << p.vy / c_cgs << ' ' << p.vz / c_cgs << ' ' << p.nu << ' ' << mpf << ' ' << p.weight << '\n';
                break;
            }
        }
    }
}

int main(int argc, char** argv) {
    std::vector<Photon> photons(0);
    size_t num_photon = load_photons("11hz/inits11.txt", photons);
    double scale_height = 4.1e15;
    double n0 = 1e10;
    double cross_section = 6.65e-25;
    double t_max = 13.7e9 * 3.1e7;
    Disk disk{10 * scale_height, scale_height, n0};

    size_t num_job = 12;

    // size_t photons_needed = 1000000;

    std::vector<std::thread> threads;

    /* size_t i = 0;
     for (; i < num_job; ++i) {
         threads.emplace_back(std::thread(evolves, i, std::ref(photons), disk, cross_section, photons_needed / num_job,
                                          t_max, std::stoi(argv[1]), std::stoi(argv[2])));
     }
 */
    size_t i = 0;
    size_t stride = num_photon / num_job;
    for (; i < num_job; ++i) {
        threads.emplace_back(
            std::thread(full_evolves, i, std::ref(photons), disk, cross_section, i * stride, (i + 1) * stride, t_max));
    }
    threads.emplace_back(
        std::thread(full_evolves, i, std::ref(photons), disk, cross_section, (i + 1) * stride, num_photon, t_max));
    for (auto& th : threads) {
        th.join();
    }
    return 0;
}