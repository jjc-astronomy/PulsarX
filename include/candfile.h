/**
 * @desc Candidate file reader shared by psrfold_fil, psrfold_fil2 and psrfold_tim
 *
 * A candidate line is a list of numbers separated by spaces or tabs. Lines
 * that do not start with a digit are not candidates.
 *
 * The meaning of the columns is taken from a header line when there is one,
 * and from the number of columns otherwise.
 *
 * 1) Header: a comment line that starts with "#id" and names every column,
 *    for example
 *
 *        #id DM accel F0 F1 F2 F3 F4 S/N
 *
 *    Names are case insensitive:
 *        id                  candidate id (not used)
 *        DM                  dispersion measure (pc/cc)
 *        F0, F1, ..., Fn     spin frequency (Hz) and its derivatives
 *                            (Hz/s, Hz/s^2, ...) at the reference epoch, to
 *                            any order: phase = sum_k Fk * t^(k+1) / (k+1)!
 *        acc | accel         line-of-sight acceleration (m/s^2)
 *        jerk                its first derivative (m/s^3)
 *        snap                its second derivative (m/s^4)
 *        crackle             its third derivative (m/s^5)
 *                            These kinematic columns give the same polynomial
 *                            as F1 ... F4: with D(t) = acc*t^2/2 + jerk*t^3/6
 *                            + snap*t^4/24 + crackle*t^5/120 the delay in
 *                            metres, phase = F0 * (t - D(t)/c), that is
 *                            F1 = -acc*F0/c, F2 = -jerk*F0/c, F3 = -snap*F0/c,
 *                            F4 = -crackle*F0/c. A kinematic value is used
 *                            when the matching Fk is not given or is zero
 *                            (as acc has always been for F1).
 *        Pb A1 T0 OM ECC     Keplerian orbit (day, lt-sec, MJD, rad, -);
 *                            with these the candidate is folded as a binary
 *                            and F0 is the intrinsic spin frequency
 *        S/N | SNR           signal-to-noise ratio of the search
 *    Any other name is a column that is skipped. A header needs at least DM
 *    and F0; one that has neither (the #id line of a psrfold output table,
 *    for instance) is not a column header and is ignored. Once a header is
 *    read, every candidate line must have as many columns as it names.
 *
 * 2) No header: by the number of columns N
 *        6          id DM acc F0 F1 S/N
 *        7          id DM acc F0 F1 F2 S/N
 *        12         id DM acc F0 F1 F2 Pb A1 T0 OM ECC S/N   (binary)
 *        8-11, >12  id DM acc F0 F1 ... F(N-5) S/N
 *    A polynomial up to F7 would also have 12 columns; it needs a header.
 */

#ifndef CANDFILE_H
#define CANDFILE_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "archivelite.h"
#include "constants.h"

namespace Pulsar
{
	class CandFile
	{
	public:
		enum ColumnType {SKIP, DM, ACC, KIN, FREQ, PB, A1, T0, OM, ECC, SNR};

		struct Column
		{
			ColumnType type;
			int order; /* k of Fk, for FREQ; k of the Fk it stands for, for KIN (jerk 2, snap 3, crackle 4) */
		};

	public:
		CandFile() : has_header(false), lineno(0) {}

		/**
		 * @brief read all candidates of a file
		 *
		 * @param filename candidate file
		 * @param base     values of the parameters a line does not give (the command line defaults)
		 * @return one ArchiveLite per candidate line, in file order
		 */
		std::vector<ArchiveLite> read(const std::string &filename, const ArchiveLite &base)
		{
			std::ifstream infile(filename);
			if (!infile.is_open())
				throw std::runtime_error("can not open candidate file " + filename);

			std::vector<ArchiveLite> cands;
			std::string line;
			lineno = 0;
			while (std::getline(infile, line))
			{
				lineno++;
				boost::trim(line);
				if (line.empty()) continue;

				if (line[0] == '#')
				{
					read_header(line);
					continue;
				}

				if (!isdigit((unsigned char)line[0])) continue;

				std::vector<std::string> items;
				boost::split(items, line, boost::is_any_of("\t "), boost::token_compress_on);

				cands.push_back(parse(items, base, filename));
			}

			return cands;
		}

		/**
		 * @brief read a header line; returns true if it defines the columns
		 */
		bool read_header(const std::string &line)
		{
			std::string s = line;
			boost::trim_left_if(s, boost::is_any_of("# \t"));

			std::vector<std::string> names;
			boost::split(names, s, boost::is_any_of("\t "), boost::token_compress_on);

			if (names.empty() || boost::to_lower_copy(names[0]) != "id") return false;

			std::vector<Column> cols;
			std::map<int, int> nfreq;
			int ndm = 0;
			for (auto name=names.begin(); name!=names.end(); ++name)
			{
				Column col = get_column(*name);
				if (col.type == DM) ndm++;
				if (col.type == FREQ) nfreq[col.order]++;
				cols.push_back(col);
			}

			if (ndm != 1 || nfreq.count(0) == 0) return false;

			for (auto f=nfreq.begin(); f!=nfreq.end(); ++f)
			{
				if (f->second != 1) return false;
			}

			layout = cols;
			has_header = true;

			return true;
		}

		bool has_column(ColumnType type) const
		{
			for (auto col=layout.begin(); col!=layout.end(); ++col)
			{
				if (col->type == type) return true;
			}
			return false;
		}

	private:
		static Column get_column(const std::string &name)
		{
			std::string s = boost::to_lower_copy(name);

			if (s == "dm") return Column{DM, 0};
			if (s == "acc" || s == "accel" || s == "acceleration") return Column{ACC, 0};
			if (s == "jerk") return Column{KIN, 2};
			if (s == "snap") return Column{KIN, 3};
			if (s == "crackle") return Column{KIN, 4};
			if (s == "pb") return Column{PB, 0};
			if (s == "a1") return Column{A1, 0};
			if (s == "t0") return Column{T0, 0};
			if (s == "om") return Column{OM, 0};
			if (s == "ecc") return Column{ECC, 0};
			if (s == "s/n" || s == "snr" || s == "sn") return Column{SNR, 0};

			if (s.size() > 1 && s[0] == 'f' && std::all_of(s.begin()+1, s.end(), [](unsigned char c){return isdigit(c);}))
				return Column{FREQ, std::stoi(s.substr(1))};

			return Column{SKIP, 0};
		}

		/**
		 * @brief column layout of a line without header, from its number of columns
		 */
		static std::vector<Column> get_layout(size_t ncol)
		{
			std::vector<Column> cols;

			if (ncol < 6) return cols;

			cols.push_back(Column{SKIP, 0});
			cols.push_back(Column{DM, 0});
			cols.push_back(Column{ACC, 0});
			cols.push_back(Column{FREQ, 0});

			if (ncol == 12)
			{
				cols.push_back(Column{FREQ, 1});
				cols.push_back(Column{FREQ, 2});
				cols.push_back(Column{PB, 0});
				cols.push_back(Column{A1, 0});
				cols.push_back(Column{T0, 0});
				cols.push_back(Column{OM, 0});
				cols.push_back(Column{ECC, 0});
			}
			else
			{
				for (size_t k=1; k+5<=ncol; k++)
					cols.push_back(Column{FREQ, (int)k});
			}

			cols.push_back(Column{SNR, 0});

			return cols;
		}

		ArchiveLite parse(const std::vector<std::string> &items, const ArchiveLite &base, const std::string &filename)
		{
			std::vector<Column> cols = has_header ? layout : get_layout(items.size());

			if (cols.size() != items.size())
			{
				std::string msg = filename + " line " + std::to_string(lineno) + ": " + std::to_string(items.size()) + " columns, ";
				if (has_header)
					msg += "the header names " + std::to_string(cols.size());
				else
					msg += "at least 6 are needed (id DM acc F0 F1 S/N)";
				throw std::runtime_error(msg);
			}

			ArchiveLite fdr = base;

			std::map<int, double> freqs;
			std::map<int, double> kins;
			double pb = 0., a1 = 0., t0 = 0., om = 0., ecc = 0.;
			bool has_pb = false, has_a1 = false, has_t0 = false;

			for (size_t k=0; k<cols.size(); k++)
			{
				if (cols[k].type == SKIP) continue;

				double value = 0.;
				try
				{
					value = std::stod(items[k]);
				}
				catch (const std::exception &e)
				{
					throw std::runtime_error(filename + " line " + std::to_string(lineno) + ": column " + std::to_string(k+1) + " (" + items[k] + ") is not a number");
				}

				switch (cols[k].type)
				{
				case DM: fdr.dm = value; break;
				case ACC: fdr.acc = value; break;
				case FREQ: freqs[cols[k].order] = value; break;
				case KIN: kins[cols[k].order] = value; break;
				case PB: pb = value; has_pb = true; break;
				case A1: a1 = value; has_a1 = true; break;
				case T0: t0 = value; has_t0 = true; break;
				case OM: om = value; break;
				case ECC: ecc = value; break;
				case SNR: fdr.snr = value; break;
				default: break;
				}
			}

			/* kinematic terms above the acceleration: Fk = -value * F0 / c, where Fk is not given or zero
			 * (the acceleration itself is applied to F1 by ArchiveLite::prepare) */
			for (auto kin=kins.begin(); kin!=kins.end(); ++kin)
			{
				if (freqs.count(kin->first) == 0 || freqs[kin->first] == 0.)
					freqs[kin->first] = -kin->second * freqs[0] / CONST_C;
			}

			/* the spin polynomial: a derivative the line does not give is zero once the line gives a higher one */
			int maxorder = freqs.rbegin()->first;
			fdr.f0 = freqs[0];
			if (maxorder >= 1) fdr.f1 = freqs.count(1) ? freqs[1] : 0.;
			if (maxorder >= 2) fdr.f2 = freqs.count(2) ? freqs[2] : 0.;
			if (maxorder >= 3)
			{
				fdr.fhigh.assign(maxorder-2, 0.);
				for (int k=3; k<=maxorder; k++)
				{
					if (freqs.count(k)) fdr.fhigh[k-3] = freqs[k];
				}
			}

			if (has_pb || has_a1 || has_t0)
			{
				if (!(has_pb && has_a1 && has_t0))
					throw std::runtime_error(filename + " line " + std::to_string(lineno) + ": a binary candidate needs Pb, A1 and T0");

				std::vector<double> kepler_params{fdr.f0, pb, a1, t0, om, ecc};
				fdr.orb = Kepler(kepler_params);
				fdr.use_kepler = true;
			}

			return fdr;
		}

	private:
		std::vector<Column> layout;
		bool has_header;
		long int lineno;
	};
}

#endif /* CANDFILE_H */
