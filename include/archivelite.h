/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 15:54:10
 * @modify date 2020-06-07 15:54:10
 * @desc [description]
 */

#ifndef ARCHIVELITE_H
#define ARCHIVELITE_H

#include <fstream>
#include <vector>
#include <string>

#include "dedispersionlite.h"
#include "predictor.h"
#include "kepler.h"
#include "constants.h"
#include "mjd.h"

using namespace std;

namespace Pulsar
{
	class IntegrationLite
	{
	public:
		IntegrationLite()
		{
			ffold = 0.;
			tsubint = 0.;
			offs_sub = 0.;
			npol = 0;
			nchan = 0;
			nbin = 0;

			mean_var_ready = false;
		}
		~IntegrationLite(){}
		void resize(int np, int nc, int nb)
		{
			npol = np;
			nchan = nc;
			nbin = nb;
			data.clear();
			data.resize(npol*nchan*nbin, 0.);
		}
	public:
		double ffold;
		double tsubint;
		double offs_sub;

		int npol;
		int nchan;
		int nbin;

		vector<float> data;
		std::vector<double> means;
		std::vector<double> vars;
		bool mean_var_ready;
	};

	class ArchiveLite
	{
	public:
		ArchiveLite();
		~ArchiveLite();
		void close()
		{
			for (auto p=profiles.begin(); p!=profiles.end(); ++p)
			{
				p->data.clear();
				p->data.shrink_to_fit();
			}
		}
		void prepare(DataBuffer<float> &databuffer);
		bool runDspsr(DataBuffer<float> &databuffer);
		bool runTRLSM(DataBuffer<float> &databuffer);
		bool runPresto(DataBuffer<float> &databuffer);
		void resize(int np, int nc, int nb);
		void dump2bin(const string &rootname)
		{
			ofstream outfile;
			outfile.open(rootname+".spft", ofstream::binary);

			for (auto pro=profiles.begin(); pro!=profiles.end(); ++pro)
			{
				outfile.write((char *)(&((*pro).data[0])), sizeof(float)*npol*nchan*nbin);
			}
			outfile.close();
		}
		void read_archive(const std::string &fname);
		/**
		 * @brief phase (turns) of the spin polynomial at t seconds from the reference epoch
		 * phase = sum_k Fk * t^(k+1) / (k+1)!, with F0, F1, F2 and the higher derivatives in fhigh
		 */
		long double get_poly_phase(long double t) const
		{
			long double phi = 0.;
			for (long int k=fhigh.size()-1; k>=0; k--)
			{
				phi = (phi + fhigh[k]) * t / (k + 4);
			}
			phi = (phi + f2) * t / 3.;
			phi = (phi + f1) * t / 2.;
			phi = (phi + f0) * t;
			return phi;
		}
		/**
		 * @brief spin frequency (Hz) of the polynomial at t seconds from the reference epoch
		 */
		long double get_poly_freq(long double t) const
		{
			long double f = 0.;
			for (long int k=fhigh.size()-1; k>=0; k--)
			{
				f = (f + fhigh[k]) * t / (k + 3);
			}
			f = (f + f2) * t / 2.;
			f = (f + f1) * t;
			f += f0;
			return f;
		}
		/**
		 * @brief Fk of the spin polynomial (F0, F1, F2 and the higher derivatives), zero beyond the last one
		 */
		double get_fn(int k) const
		{
			if (k == 0) return f0;
			if (k == 1) return f1;
			if (k == 2) return f2;
			if (k >= 3 && k-3 < (long int)fhigh.size()) return fhigh[k-3];
			return 0.;
		}
		/**
		 * @brief order of the highest non-zero derivative of the spin polynomial (0 for a constant frequency)
		 */
		int get_max_order() const
		{
			for (int k=fhigh.size()+2; k>=1; k--)
			{
				if (get_fn(k) != 0.) return k;
			}
			return 0;
		}
		/**
		 * @brief the kinematic equivalent of Fk: the (k-1)th time derivative of the line-of-sight
		 * acceleration in m/s^(k+1), -Fk/F0*c (k = 1 acceleration, 2 jerk, 3 snap, 4 crackle)
		 */
		double get_kinematic(int k) const
		{
			if (f0 == 0.) return 0.;
			return -get_fn(k) / f0 * CONST_C;
		}
		static std::string get_kinematic_name(int k)
		{
			switch (k)
			{
			case 1: return "acc";
			case 2: return "jerk";
			case 3: return "snap";
			case 4: return "crackle";
			case 5: return "pop";
			default: return "kin" + std::to_string(k);
			}
		}
	private:
		double get_phase(MJD mjd, MJD &mjdref)
		{
			if (use_t2pred)
			{
				return pred.get_fphase(mjd.to_day(), fref);
			}
			else if (use_kepler)
			{
				return orb.get_fphase(mjd.to_day(), mjdref.to_day());
			}
			else
			{
				long double t = (mjd-mjdref).to_second();
				double phi = get_poly_phase(t);
				phi -= floor(phi);
				return phi;
			}
		}
		double get_ffold(MJD mjd, MJD &mjdref)
		{
			if (use_t2pred)
			{
				double pfold = pred.get_pfold(mjd.to_day(), fref);
				return 1./pfold;
			}
			else if (use_kepler)
			{
				return orb.get_ffold(mjd.to_day());
			}
			else
			{
				long double t = (mjd-mjdref).to_second();
				double f = get_poly_freq(t);
				return f;
			}
		}
		MJD get_epoch(MJD &start_time, MJD &end_time, MJD &ref_epoch)
		{
			long double t = (start_time.to_second()+end_time.to_second())/2.;
			t -= ref_epoch.to_second();
			long int phi = floor(get_poly_phase(t));
			long double tepoch = 0.;
			if (!fhigh.empty())
			{
				/* time near t where the polynomial phase is the integer phi (Newton) */
				int n = 16;
				tepoch = t;
				long double dt = 1.;
				while (std::abs(dt) > 1e-9 && --n >= 0)
				{
					long double f = get_poly_freq(tepoch);
					if (f == 0.) break;
					dt = (phi - get_poly_phase(tepoch)) / f;
					tepoch += dt;
				}
			}
			else if (f1 != 0 || f2 != 0)
			{
				if (f2 == 0)
				{
					double tmp = f0*f0+2*f1*phi;
					tmp = tmp<0? 0:tmp;
					long double tepoch1 = (sqrt(tmp)-f0)/f1;
					long double tepoch2 = (-sqrt(f0*f0+2*f1*phi)-f0)/f1;
					tepoch = abs(t-tepoch1)<abs(t-tepoch2) ? tepoch1:tepoch2;
				}
				else
				{
					double a = phi / f0;
					double b = -0.5 * f1 / f0;
					double c = -f2 / f0 / 6.;

					int n = 8;
					tepoch = t;
					long double dt = 1.;
					while (std::abs(dt) > 1e-9 && --n >= 0)
					{
						dt = a - tepoch  + b * tepoch * tepoch + c * tepoch * tepoch * tepoch;
						tepoch += dt;
					}
				}
			}
			else
				tepoch = phi/f0;

			return MJD((tepoch+ref_epoch.to_second())/86400.);
		}
	public:
		MJD start_mjd;
		MJD ref_epoch;
		double f0;
		double f1;
		double f2;
		std::vector<double> fhigh; /* F3, F4, ... (Hz/s^3, Hz/s^4, ...) */
		double acc;
		double dm;
		double snr;
		long int nblock;
		bool use_t2pred;
		Predictors pred;
		bool use_kepler;
		Kepler orb;
		bool dedispersed;
	public:
		int nbin;
		int nchan;
		int npol;
		double tbin;
		vector<double> frequencies;
		vector<IntegrationLite> profiles;
		MJD sub_mjd;
		IntegrationLite sub_int;
		std::vector<double> means;
		std::vector<double> vars;
		bool mean_var_ready;
	private:
		double fref;
		int iblock;
		vector<float> mxWTW;
		vector<float> vWTd_T;
		vector<int> hits;
		vector<float> profilesTPF;
		vector<float> profilesPFT;
	};
}

#endif /* ARCHIVELITE_H */
