#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

struct SNP { long long pos1; std::string id, ref, alt; double af; };
struct LDPair { int a,b,bin; };
static std::vector<std::string> split(const std::string& s, char d='\t'){ std::vector<std::string> v; std::string x; std::stringstream ss(s); while(std::getline(ss,x,d)) v.push_back(x); return v; }
static uint64_t mix64(uint64_t x){ x += 0x9e3779b97f4a7c15ULL; x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL; x=(x^(x>>27))*0x94d049bb133111ebULL; return x^(x>>31); }
static uint64_t keyseed(uint64_t base, uint64_t a, uint64_t b=0){ return mix64(base ^ mix64(a+0x123456789ULL) ^ mix64(b+0xABCDEFULL)); }
static std::vector<SNP> read_snps(const std::string& path){
    std::ifstream in(path); if(!in){ std::cerr<<"Cannot open "<<path<<"\n"; std::exit(2);} std::string line; std::getline(in,line); std::vector<SNP> s;
    while(std::getline(in,line)){ if(line.empty()) continue; auto f=split(line); if(f.size()<10) continue; s.push_back({std::stoll(f[2]),f[4],f[5],f[6],std::stod(f[9])}); }
    return s;
}
static std::vector<LDPair> read_ldpairs(const std::string& path){
    std::ifstream in(path); if(!in) return {}; std::string line; std::getline(in,line); std::vector<LDPair> v;
    while(std::getline(in,line)){ if(line.empty()) continue; auto f=split(line); if(f.size()<4) continue; v.push_back({std::stoi(f[0]),std::stoi(f[1]),std::stoi(f[3])}); } return v;
}
static std::vector<int> parse_int_csv(const std::string& s){ std::vector<int> v; std::stringstream ss(s); std::string x; while(std::getline(ss,x,',')) if(!x.empty()) v.push_back(std::stoi(x)); std::sort(v.begin(),v.end()); return v; }

struct RecombSegment { long long first_threshold1, last_threshold1; };
struct RecombMap {
    // Poisson mean integrates rate over all [0, chromosome_length) bases.
    double expected_crossovers=0.0;
    double placement_mass=0.0;
    std::vector<RecombSegment> segments;
    std::vector<double> cumulative_mass;
};

static RecombMap read_recomb_map(const std::string& path, long long chr_len){
    if(chr_len<2) throw std::runtime_error("Recombination map requires chromosome length >= 2");
    std::ifstream in(path);
    if(!in) throw std::runtime_error("Cannot open recombination map: "+path);
    std::string line;
    if(!std::getline(in,line)) throw std::runtime_error("Empty recombination map: "+path);
    auto header=split(line);
    if(header.size()<3 || header[0]!="start0" || header[1]!="end0" || header[2]!="rate_per_bp")
        throw std::runtime_error("Recombination map header must begin start0,end0,rate_per_bp (tab separated)");
    RecombMap result;
    long long previous_end=0;
    size_t row=1;
    while(std::getline(in,line)){
        ++row;
        if(line.empty()) continue;
        auto f=split(line);
        if(f.size()<3) throw std::runtime_error("Incomplete recombination map row "+std::to_string(row));
        size_t ns=0,ne=0,nr=0;
        long long start=std::stoll(f[0],&ns),end=std::stoll(f[1],&ne);
        double rate=std::stod(f[2],&nr);
        if(ns!=f[0].size() || ne!=f[1].size() || nr!=f[2].size())
            throw std::runtime_error("Malformed number in recombination map row "+std::to_string(row));
        if(start!=previous_end || end<=start || end>chr_len)
            throw std::runtime_error("Recombination map must be ordered, contiguous and cover exactly [0, chrlen); row "+std::to_string(row));
        if(!std::isfinite(rate) || rate<0)
            throw std::runtime_error("Recombination rate must be finite and nonnegative; row "+std::to_string(row));
        result.expected_crossovers+=(double)(end-start)*rate;
        // Preserve the legacy crossover threshold x in [1, chrlen-1].
        // SNPs with POS < x come from the current homolog, POS >= x from the
        // next homolog. The rate at BED coordinate x-1 weights this threshold.
        // Thus the last terminal base is not an admissible threshold. Event
        // counts still integrate the complete chromosome (legacy rate * L);
        // placement is conditional on the admissible thresholds. No clamping
        // or extrapolation of sampled crossovers at chromosome ends is used.
        long long clipped_end=std::min(end,chr_len-1);
        if(rate>0 && clipped_end>start){
            double mass=(double)(clipped_end-start)*rate;
            double next_mass=result.placement_mass+mass;
            if(!std::isfinite(next_mass) || next_mass<=result.placement_mass)
                throw std::runtime_error("Recombination placement weights overflow or lose numerical precision");
            result.placement_mass=next_mass;
            result.segments.push_back({start+1,clipped_end});
            result.cumulative_mass.push_back(result.placement_mass);
        }
        previous_end=end;
    }
    if(previous_end!=chr_len)
        throw std::runtime_error("Recombination map does not cover the complete chromosome");
    if(!std::isfinite(result.expected_crossovers) || result.expected_crossovers<=0 || result.placement_mass<=0)
        throw std::runtime_error("Recombination map needs positive, finite total mass and admissible placement mass");
    return result;
}

static long long sample_recomb_threshold(const RecombMap& map,std::mt19937_64& rng){
    std::uniform_real_distribution<double> draw(0.0,map.placement_mass);
    double target=std::min(draw(rng),std::nextafter(map.placement_mass,0.0));
    size_t i=(size_t)(std::upper_bound(map.cumulative_mass.begin(),map.cumulative_mass.end(),target)-map.cumulative_mass.begin());
    const auto& segment=map.segments.at(i);
    std::uniform_int_distribution<long long> position(segment.first_threshold1,segment.last_threshold1);
    return position(rng);
}

struct Pop {
    int N; size_t S; std::vector<uint8_t> h;
    Pop(int n,size_t s):N(n),S(s),h((size_t)2*n*s,0){}
    uint8_t* row(int r){ return h.data()+(size_t)r*S; }
    const uint8_t* row(int r) const { return h.data()+(size_t)r*S; }
};

// LD-free founders: each SNP preserves round(AF * 2N) ALT copies exactly, but copies are independently shuffled across haplotypes.
// RNG is keyed by SNP index, so output is invariant to OpenMP thread count.
static void init_founders(Pop& pop, const std::vector<SNP>& snps, uint64_t seed, int threads){
    const int H=2*pop.N; const size_t S=pop.S;
#ifdef _OPENMP
    omp_set_num_threads(threads);
#pragma omp parallel for schedule(static)
#endif
    for(long long j=0;j<(long long)S;j++){
        std::mt19937_64 rng(keyseed(seed,0xF00DULL,(uint64_t)j));
        std::vector<int> idx(H); std::iota(idx.begin(),idx.end(),0); std::shuffle(idx.begin(),idx.end(),rng);
        int k=(int)std::llround(snps[(size_t)j].af*H); if(k<0)k=0; if(k>H)k=H;
        for(int z=0;z<k;z++) pop.h[(size_t)idx[z]*S+(size_t)j]=1;
    }
}

static void make_gamete(const Pop& pop, int parent, uint8_t* dst, const std::vector<long long>& pos, double recomb_rate, long long chr_len, std::mt19937_64& rng, const RecombMap* recomb_map=nullptr){
    const size_t S=pop.S; std::bernoulli_distribution coin(0.5); int hap=coin(rng)?1:0;
    double lambda=recomb_map?recomb_map->expected_crossovers:recomb_rate*(double)chr_len; std::poisson_distribution<int> P(lambda); int nxo=P(rng);
    std::vector<long long> xo; xo.reserve(nxo);
    if(nxo>0){
        if(recomb_map){ for(int z=0;z<nxo;z++) xo.push_back(sample_recomb_threshold(*recomb_map,rng)); }
        else { std::uniform_int_distribution<long long> X(1,std::max(1LL,chr_len-1)); for(int z=0;z<nxo;z++) xo.push_back(X(rng)); }
        std::sort(xo.begin(),xo.end());
    }
    size_t start=0;
    for(long long x: xo){ size_t end=std::lower_bound(pos.begin()+start,pos.end(),x)-pos.begin(); if(end>start){ const uint8_t* src=pop.row(2*parent+hap); std::copy(src+start,src+end,dst+start);} start=end; hap^=1; }
    if(start<S){ const uint8_t* src=pop.row(2*parent+hap); std::copy(src+start,src+S,dst+start); }
}

// Pure neutral Wright-Fisher-style reproduction. No agLink input or fitness term exists.
// RNG is keyed by generation and offspring index, so output is invariant to task scheduling/thread count.
static void next_generation(const Pop& pop, Pop& nxt, const std::vector<long long>& pos, double selfing, double recomb_rate, long long chr_len, uint64_t seed, int gen, int threads, const RecombMap* recomb_map=nullptr){
#ifdef _OPENMP
    omp_set_num_threads(threads);
#pragma omp parallel for schedule(static)
#endif
    for(int i=0;i<pop.N;i++){
        std::mt19937_64 rng(keyseed(seed,(uint64_t)gen,(uint64_t)i));
        std::uniform_int_distribution<int> PARENT(0,pop.N-1); std::bernoulli_distribution SELF(selfing);
        int p1=PARENT(rng); int p2=SELF(rng)?p1:PARENT(rng);
        make_gamete(pop,p1,nxt.row(2*i),pos,recomb_rate,chr_len,rng,recomb_map);
        make_gamete(pop,p2,nxt.row(2*i+1),pos,recomb_rate,chr_len,rng,recomb_map);
    }
}

static std::vector<int> sample_indices(const Pop& pop, int outN, uint64_t seed, int gen){
    outN=std::min(outN,pop.N); std::vector<int> inds(pop.N); std::iota(inds.begin(),inds.end(),0);
    std::mt19937_64 rng(keyseed(seed,0x5A5AULL,(uint64_t)gen)); std::shuffle(inds.begin(),inds.end(),rng); inds.resize(outN); std::sort(inds.begin(),inds.end()); return inds;
}

static void write_qc(std::ofstream& ldout, std::ofstream& mafout, int gen, const Pop& pop, const std::vector<SNP>& snps, const std::vector<LDPair>& pairs, int outN, uint64_t seed){
    auto inds=sample_indices(pop,outN,seed,gen);
    int maxbin=-1; for(const auto&p:pairs) maxbin=std::max(maxbin,p.bin); std::vector<double> sum(maxbin+1,0.0); std::vector<long long> cnt(maxbin+1,0);
    for(const auto&p:pairs){
        double sx=0,sy=0,sxx=0,syy=0,sxy=0; int n=0;
        for(int i:inds){ double x=(double)pop.row(2*i)[p.a]+pop.row(2*i+1)[p.a]; double y=(double)pop.row(2*i)[p.b]+pop.row(2*i+1)[p.b]; sx+=x;sy+=y;sxx+=x*x;syy+=y*y;sxy+=x*y;n++; }
        if(n<3) continue; double vx=sxx-sx*sx/n, vy=syy-sy*sy/n; if(vx<=0||vy<=0) continue; double cov=sxy-sx*sy/n; double r2=(cov*cov)/(vx*vy); if(std::isfinite(r2)){ sum[p.bin]+=r2; cnt[p.bin]++; }
    }
    for(int b=0;b<=maxbin;b++) ldout<<gen<<"\t"<<b<<"\t"<<cnt[b]<<"\t"<<(cnt[b]?sum[b]/cnt[b]:NAN)<<"\n";
    size_t step=std::max<size_t>(1, snps.size()/20000); double se=0,sa=0,ho=0,hs=0; long long nf=0,nm=0;
    for(size_t j=0;j<snps.size();j+=step){ long long ac=0; for(int i:inds) ac += pop.row(2*i)[j]+pop.row(2*i+1)[j]; double q=(double)ac/(2.0*inds.size()); double p=snps[j].af; double d=q-p; se+=d*d; sa+=std::abs(d); ho+=2*p*(1-p); hs+=2*q*(1-q); nf += (q<=0||q>=1); nm++; }
    mafout<<gen<<"\t"<<nm<<"\t"<<std::sqrt(se/nm)<<"\t"<<sa/nm<<"\t"<<(double)nf/nm<<"\t"<<ho/nm<<"\t"<<hs/nm<<"\n";
}

static void write_vcf(const std::string& path, const std::string& chrom, const Pop& pop, const std::vector<SNP>& snps, int outN, uint64_t seed, int gen){
    auto inds=sample_indices(pop,outN,seed,gen); std::ofstream out(path); if(!out){ std::cerr<<"Cannot write "<<path<<"\n"; std::exit(3);}
    out<<"##fileformat=VCFv4.2\n##source=aglink_neutral_forward_v19\n##contig=<ID="<<chrom<<">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n";
    out<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"; for(int i:inds) out<<"\tSIM"<<i+1; out<<"\n";
    for(size_t j=0;j<snps.size();j++){ const auto&s=snps[j]; out<<chrom<<"\t"<<s.pos1<<"\t"<<s.id<<"\t"<<s.ref<<"\t"<<s.alt<<"\t.\tPASS\t.\tGT"; for(int i:inds) out<<"\t"<<(int)pop.row(2*i)[j]<<"|"<<(int)pop.row(2*i+1)[j]; out<<"\n"; }
}

int main(int argc,char**argv){
    std::string snpfile,outdir,chrom="1A",checkcsv="10,20,50,100,200",ldpairfile,ldoutpath,mafoutpath,recombmapfile; int N=203,outN=203,threads=1,emit_vcf=1; double recomb=1e-8,selfing=0.95; long long chrlen=0; uint64_t seed=1;
    for(int i=1;i<argc;i++){
        std::string a=argv[i];
        auto need=[&](){ if(i+1>=argc){std::cerr<<"Missing value for "<<a<<"\n";std::exit(1);} return std::string(argv[++i]);};
        if(a=="--snps")snpfile=need();
        else if(a=="--outdir")outdir=need();
        else if(a=="--chrom")chrom=need();
        else if(a=="--N")N=std::stoi(need());
        else if(a=="--outN")outN=std::stoi(need());
        else if(a=="--recomb")recomb=std::stod(need());
        else if(a=="--recomb-map")recombmapfile=need();
        else if(a=="--selfing")selfing=std::stod(need());
        else if(a=="--chrlen")chrlen=std::stoll(need());
        else if(a=="--seed")seed=std::stoull(need());
        else if(a=="--threads")threads=std::stoi(need());
        else if(a=="--checkpoints")checkcsv=need();
        else if(a=="--emit-vcf")emit_vcf=std::stoi(need());
        else if(a=="--ldpairs")ldpairfile=need();
        else if(a=="--ld-out")ldoutpath=need();
        else if(a=="--maf-out")mafoutpath=need();
    }
    if(snpfile.empty()||outdir.empty()){ std::cerr<<"Required: --snps --outdir\n"; return 1; }
    auto snps=read_snps(snpfile); if(snps.empty()){ std::cerr<<"No SNPs\n"; return 2; } if(chrlen<=0) chrlen=snps.back().pos1;
    RecombMap spatial_map;
    const RecombMap* spatial_map_ptr=nullptr;
    if(!recombmapfile.empty()){
        try {
            spatial_map=read_recomb_map(recombmapfile,chrlen);
            for(const auto& snp:snps) if(snp.pos1<1 || snp.pos1>chrlen)
                throw std::runtime_error("SNP coordinate lies outside the recombination map chromosome");
        } catch(const std::exception& error){
            std::cerr<<"Invalid recombination map: "<<error.what()<<"\n";
            return 4;
        }
        spatial_map_ptr=&spatial_map;
        std::cerr<<"Spatial recombination map="<<recombmapfile
                 <<" expected_crossovers_per_gamete="<<spatial_map.expected_crossovers
                 <<" mean_rate_per_bp="<<spatial_map.expected_crossovers/(double)chrlen
                 <<" positive_segments="<<spatial_map.segments.size()<<"\n";
    }
    auto ldpairs=ldpairfile.empty()?std::vector<LDPair>{}:read_ldpairs(ldpairfile); for(auto&p:ldpairs) if(p.a<0||p.b<0||p.a>=(int)snps.size()||p.b>=(int)snps.size()){std::cerr<<"Bad LD pair index\n";return 5;}
    std::ofstream ldout,mafout; if(!ldoutpath.empty()){ldout.open(ldoutpath);ldout<<"generation\tbin_index\tn_pairs\tmean_r2\n";} if(!mafoutpath.empty()){mafout.open(mafoutpath);mafout<<"generation\tn_snps_sampled\taf_rmse\taf_mean_abs_diff\tfixed_fraction\tmean_het_G0_target\tmean_het_simulated\n";}
    std::vector<long long> pos; pos.reserve(snps.size()); for(auto&s:snps)pos.push_back(s.pos1); auto cps=parse_int_csv(checkcsv); if(cps.empty()) return 6; int maxgen=cps.back();
    Pop pop(N,snps.size()),nxt(N,snps.size()); init_founders(pop,snps,seed,threads);
    auto emit=[&](int g){ if(ldout.is_open()&&mafout.is_open()) write_qc(ldout,mafout,g,pop,snps,ldpairs,outN,seed); if(emit_vcf) write_vcf(outdir+"/G"+std::to_string(g)+".vcf",chrom,pop,snps,outN,seed,g); };
    if(std::binary_search(cps.begin(),cps.end(),0)) emit(0);
    for(int g=1;g<=maxgen;g++){ next_generation(pop,nxt,pos,selfing,recomb,chrlen,seed,g,threads,spatial_map_ptr); pop.h.swap(nxt.h); if(std::binary_search(cps.begin(),cps.end(),g)) emit(g); }
    std::cerr<<"Done neutral chrom="<<chrom<<" N="<<N<<" SNPs="<<snps.size()<<" maxgen="<<maxgen<<"\n"; return 0;
}
