#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include "cmdline.h"
#include "vcf.h"


using std::size_t;


// TODO
// Chi-squared test for independance
// void calc_chisq(const std::vector<allele_t> &x, const std::vector<allele_t> &y);


namespace {


struct Parameter
{
    std::string vcf;
    std::string out;
    std::string loc;
    double rsq = 0.5;
    int maxdist = 500000;
} par ;


// LD (D')  -  Lewontin, R.C. (1964). Genetics 49(1), 49-67.
// LD (r^2) -  Hill, W.G., and Robertson, A. (1968). Theor Appl Genet 38(6), 226-231.

void calc_dprime_rsq_kernel(double pa, double pb, double pab, double &dprime, double &rsq)
{
    auto D = pab - pa*pb;
    auto Dmax = D > 0 ? std::min(pa*(1-pb),(1-pa)*pb) : std::min(pa*pb,(1-pa)*(1-pb));
    dprime = std::fabs(D) / Dmax;  // fabs is required for multi-allelic D'
    rsq = D*D / (pa * (1-pa) * pb * (1-pb));
}

void calc_dprime_rsq_hap2(const std::vector<allele_t> &x, const std::vector<allele_t> &y, double &dprime, double &rsq)
{
    static const allele_t a = 1, b = 1;

    auto n = x.size();
    int nt = 0, na = 0, nb = 0, nab = 0;

    for (size_t i = 0; i < n; ++i) {
        if (x[i] && y[i]) {
            ++nt;
            if (x[i] == a) {
                ++na;
                if (y[i] == b)
                    ++nab;
            }
            if (y[i] == b)
                ++nb;
        }
    }

    dprime = rsq = std::numeric_limits<double>::quiet_NaN();

    if (nt == 0)
        return;

    double fnt = nt;

    calc_dprime_rsq_kernel(na/fnt, nb/fnt, nab/fnt, dprime, rsq);
}

void calc_dprime_rsq_hap(const std::vector<allele_t> &x, const std::vector<allele_t> &y, double &dprime, double &rsq)
{
    auto n = x.size();

    auto xmax = * std::max_element(x.begin(), x.end());
    auto ymax = * std::max_element(y.begin(), y.end());

    dprime = rsq = std::numeric_limits<double>::quiet_NaN();

    if (xmax < 2 || ymax < 2)
        return;

    if (xmax == 2 && ymax == 2) {
        calc_dprime_rsq_hap2(x, y, dprime, rsq);
        return;
    }

    std::vector<size_t> idx;
    idx.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        if (x[i] && y[i])
            idx.push_back(i);
    }

    if ( idx.empty() )
        return;

    dprime = rsq = 0.0;

    for (allele_t a = 1; a <= xmax; ++a) {
        for (allele_t b = 1; b <= ymax; ++b) {
            int na = 0, nb = 0, nab = 0;

            for (auto i : idx) {
                if (x[i] == a) {
                    ++na;
                    if (y[i] == b)
                        ++nab;
                }
                if (y[i] == b)
                    ++nb;
            }

            if (na == 0 || nb == 0)
                continue;

            double nt = idx.size();
            double pa = na / nt;
            double pb = nb / nt;
            double pab = nab / nt;

            double t1 = 0.0, t2 = 0.0;
            calc_dprime_rsq_kernel(pa, pb, pab, t1, t2);

            pab = pa * pb;
            dprime += pab * t1;
            rsq += pab * t2;
        }
    }
}

// Estimates haplotype frequencies via the EM algorithm
// AB, Ab, aB, ab, AaBb
void calc_hap_prob_EM(int n11, int n12, int n21, int n22, int ndh, double &p11, double &p12, double &p21, double &p22)
{
    static const int maxit = 1000;
    static const double tol = 1e-10;

    double n = n11 + n12 + n21 + n22 + ndh * 2;
    p11 = n11 / n;
    p12 = n12 / n;
    p21 = n21 / n;
    p22 = n22 / n;

    if (ndh == 0)
        return;

    auto cp11 = p11;
    auto cp12 = p12;
    auto cp21 = p21;
    auto cp22 = p22;

    auto h = ndh / n;
    auto x = h / 2;
    auto y = h - x;

    for (int i = 0; i < maxit; ++i) {
        p11 = cp11 + x;
        p12 = cp12 + y;
        p21 = cp21 + y;
        p22 = cp22 + x;
        auto z = h * p11 * p22 / (p11 * p22 + p12 * p21);
        if (std::fabs(x - z) < tol)
            break;
        x = z;
        y = h - x;
    }
}

void calc_dprime_rsq_dip2(const std::vector<allele_t> &x, const std::vector<allele_t> &y, double &dprime, double &rsq)
{
    auto n = x.size() / 2;

    int nt = 0;
    int f[3][3] = { { 0,0,0 }, { 0,0,0 }, { 0,0,0 } };

    for (size_t i = 0; i < n; ++i) {
        auto j = i*2, k = i*2+1;
        if (x[j] && x[k] && y[j] && y[k]) {
            auto a = x[j] + x[k] - 2;
            auto b = y[j] + y[k] - 2;
            ++f[a][b];
            ++nt;
        }
    }

    dprime = rsq = std::numeric_limits<double>::quiet_NaN();

    if (nt == 0)
        return;

    int n11 = f[0][0] * 2 + f[0][1] + f[1][0];
    int n12 = f[0][2] * 2 + f[0][1] + f[1][2];
    int n21 = f[2][0] * 2 + f[1][0] + f[2][1];
    int n22 = f[2][2] * 2 + f[2][1] + f[1][2];
    int ndh = f[1][1];

    double p11 = 0.0, p12 = 0.0, p21 = 0.0, p22 = 0.0;
    calc_hap_prob_EM(n11, n12, n21, n22, ndh, p11, p12, p21, p22);

    calc_dprime_rsq_kernel(p11 + p12, p11 + p21, p11, dprime, rsq);
}

// presuming that gametic phase is known
void calc_dprime_rsq_dip(const std::vector<allele_t> &x, const std::vector<allele_t> &y, double &dprime, double &rsq)
{
    auto n = x.size() / 2;

    auto xmax = * std::max_element(x.begin(), x.end());
    auto ymax = * std::max_element(y.begin(), y.end());

    dprime = rsq = std::numeric_limits<double>::quiet_NaN();

    if (xmax < 2 || ymax < 2)
        return;

    if (xmax == 2 && ymax == 2) {
        calc_dprime_rsq_dip2(x, y, dprime, rsq);
        return;
    }

    std::vector<size_t> idx;
    idx.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        auto j = i*2, k = i*2+1;
        if (x[j] && x[k] && y[j] && y[k])
            idx.push_back(i);
    }

    if ( idx.empty() )
        return;

    dprime = rsq = 0.0;

    for (allele_t a = 1; a <= xmax; ++a) {
        for (allele_t b = 1; b <= ymax; ++b) {
            int na = 0, nb = 0, nab = 0;

            for (auto i : idx) {
                auto j = i*2, k = i*2+1;
                if (x[j] == a) {
                    ++na;
                    if (y[j] == b)
                        ++nab;
                }

                if (x[k] == a) {
                    ++na;
                    if (y[k] == b)
                        ++nab;
                }

                if (y[j] == b)
                    ++nb;

                if (y[k] == b)
                    ++nb;
            }

            if (na == 0 || nb == 0)
                continue;

            double nt = idx.size() * 2;
            auto pa = na / nt;
            auto pb = nb / nt;
            auto pab = nab / nt;

            double t1 = 0.0, t2 = 0.0;
            calc_dprime_rsq_kernel(pa, pb, pab, t1, t2);

            pab = pa * pb;
            dprime += pab * t1;
            rsq += pab * t2;
        }
    }
}

std::vector<std::string> read_string_list(const std::string &filename)
{
    std::vector<std::string> vs;

    std::ifstream ifs(filename);

    if ( ! ifs )
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
    else
        std::copy(std::istream_iterator<std::string>(ifs), std::istream_iterator<std::string>(),
                  std::back_inserter(vs));

    return vs;
}

int ldstat_list(const Genotype &gt)
{
    auto loc = read_string_list(par.loc);

    std::sort(loc.begin(), loc.end());

    std::ofstream ofs(par.out + ".sub");
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << par.out << ".list\n";
        return 1;
    }

    ofs << "Locus1\tChromosome1\tPosition1\tLocus2\tChromosome2\tPosition2\tDPrime\tRSquare\n";

    auto m = gt.loc.size();

    for (size_t i = 0; i < m; ++i) {
        if ( ! std::binary_search(loc.begin(), loc.end(), gt.loc[i]) )
            continue;

        std::cerr << "INFO: locus " << gt.loc[i] << "\n";

        for (size_t j = 0; j < m; ++j) {
            if (gt.loc[j] == gt.loc[i])
                continue;

            auto n1 = gt.allele[i].size();
            auto n2 = gt.allele[j].size();

            if (n1 < 2 || n2 < 2)
                continue;

            auto dprime = std::numeric_limits<double>::quiet_NaN();
            auto rsq = std::numeric_limits<double>::quiet_NaN();

            if (n1 == 2 && n2 == 2) {
                if (gt.ploidy == 1)
                    calc_dprime_rsq_hap2(gt.dat[i], gt.dat[j], dprime, rsq);
                else
                    calc_dprime_rsq_dip2(gt.dat[i], gt.dat[j], dprime, rsq);
            }
            else {
                if (gt.ploidy == 1)
                    calc_dprime_rsq_hap(gt.dat[i], gt.dat[j], dprime, rsq);
                else
                    calc_dprime_rsq_dip(gt.dat[i], gt.dat[j], dprime, rsq);
            }

            if ( ! std::isfinite(dprime) || ! std::isfinite(rsq) )
                continue;

            if (rsq < par.rsq)
                continue;

            ofs << gt.loc[i] << "\t" << gt.chr[i] << "\t" << gt.pos[i] << "\t"
                << gt.loc[j] << "\t" << gt.chr[j] << "\t" << gt.pos[j] << "\t"
                << dprime << "\t" << rsq << "\n";
        }
    }

    return 0;
}


} // namespace


int ldstat(int argc, char *argv[])
{
    std::cerr << "LDSTAT (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd;

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--loc", "locus list file", "");
    cmd.add("--loc-min-r2", "minimum LD (r2) threshold for locus list", "0.5");
    cmd.add("--out", "output file", "ldstat.out");
    cmd.add("--maxdist", "maximum inter-variant distance", "500000");

    cmd.parse(argc, argv);

    if (argc < 2) {
        cmd.show();
        return 1;
    }

    par.vcf = cmd.get("--vcf");
    par.loc = cmd.get("--loc");
    par.out = cmd.get("--out");
    par.rsq = std::stod(cmd.get("--loc-min-r2"));
    par.maxdist = std::stoi(cmd.get("--maxdist"));

    Genotype gt;

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 2;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    if ( ! par.loc.empty() )
        return ldstat_list(gt);

    std::ofstream ofs1(par.out + ".all");
    if ( ! ofs1 ) {
        std::cerr << "ERROR: can't open file: " << par.out << ".all\n";
        return 1;
    }

    ofs1 << "Chromosome\tLocus1\tPosition1\tLocus2\tPosition2\tDistance\tDPrime\tRSquare\n";

    int n = par.maxdist / 1000;
    using namespace boost::accumulators;
    std::vector< accumulator_set<double, stats<tag::count, tag::mean> > > acc1(n), acc2(n);

    auto m = gt.loc.size();

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            if (gt.chr[j] != gt.chr[i])
                continue;

            auto dist = std::abs(gt.pos[j] - gt.pos[i]);
            if (dist > par.maxdist)
                continue;

            auto n1 = gt.allele[i].size();
            auto n2 = gt.allele[j].size();

            if (n1 < 2 || n2 < 2)
                continue;

            auto dprime = std::numeric_limits<double>::quiet_NaN();
            auto rsq = std::numeric_limits<double>::quiet_NaN();

            if (n1 == 2 && n2 == 2) {
                if (gt.ploidy == 1)
                    calc_dprime_rsq_hap2(gt.dat[i], gt.dat[j], dprime, rsq);
                else
                    calc_dprime_rsq_dip2(gt.dat[i], gt.dat[j], dprime, rsq);
            }
            else {
                if (gt.ploidy == 1)
                    calc_dprime_rsq_hap(gt.dat[i], gt.dat[j], dprime, rsq);
                else
                    calc_dprime_rsq_dip(gt.dat[i], gt.dat[j], dprime, rsq);
            }

            if ( ! std::isfinite(dprime) || ! std::isfinite(rsq) )
                continue;

            int k = dist / 1000;
            if (dist % 1000 == 0)
                --k;

            acc1[k](dprime);
            acc2[k](rsq);

            ofs1 << gt.chr[i] << "\t"
                 << gt.loc[i] << "\t" << gt.pos[i] << "\t"
                 << gt.loc[j] << "\t" << gt.pos[j] << "\t"
                 << dist << "\t" << dprime << "\t" << rsq << "\n";
        }
    }

    std::ofstream ofs2(par.out + ".sum");
    if ( ! ofs2 ) {
        std::cerr << "ERROR: can't open file: " << par.out << ".sum\n";
        return 1;
    }

    ofs2 << "Group\tCount\tDPrime\tRSquare\n";

    for (int i = 0; i < n; ++i) {
        auto c = count(acc1[i]);
        ofs2 << (i + 1) * 1000 << "\t" << c << "\t";
        if (c != 0)
            ofs2 << mean(acc1[i]) << "\t" << mean(acc2[i]) << "\n";
        else
            ofs2 << "NA\tNA\n";
    }

    return 0;
}
