#ifndef SLICE2DFIT_C
#define SLICE2DFIT_C

#include "slice2D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>
#include <cerrno>
#include <sys/stat.h>
#include <sys/types.h>

static inline bool InsideEllipse2D(double x, double y,
                                 double xc, double yc,
                                 double rx, double ry)
{
  if (rx <= 0.0 || ry <= 0.0) return false;
  const double dx = (x - xc) / rx;
  const double dy = (y - yc) / ry;
  return (dx * dx + dy * dy) <= 1.0;
}

static inline bool EllipseOverlap2D(double x1, double y1,
                                  double x2, double y2,
                                  double rx, double ry)
{
  if (rx <= 0.0 || ry <= 0.0) return false;
  const double dx = (x1 - x2) / rx;
  const double dy = (y1 - y2) / ry;
  return (dx * dx + dy * dy) <= 4.0;
}

// Supplement the purely geometric overlap test with the experimental spectrum.
// For two nearby, non-overlapping extraction ellipses, sample DI at the point
// where the line joining the peaks crosses each ellipse.  If the observed
// intensity is still at least half the experimental intensity at the other
// peak centre, the neighbouring peak is materially contributing at the ROI
// boundary and the two fit groups should be merged.
static inline bool EdgeIntensityOverlap2D(const slice2D& s,
                                          const FitPeak2DLocal& a,
                                          const FitPeak2DLocal& b,
                                          double rx, double ry)
{
  if (rx <= 0.0 || ry <= 0.0 || s.si < 1 || s.sj < 1) return false;

  const double dx = b.x - a.x;
  const double dy = b.y - a.y;
  const double dnorm = std::sqrt((dx * dx) / (rx * rx) +
                                 (dy * dy) / (ry * ry));

  // dnorm <= 2 is already handled by EllipseOverlap2D.  Limit this
  // intensity-based extension to genuinely nearby peaks (four ROI radii).
  if (!(dnorm > 2.0) || dnorm > 4.0) return false;

  const double stepI = (s.si > 1) ? (s.ivals[1] - s.ivals[0]) : 1.0;
  const double stepJ = (s.sj > 1) ? (s.jvals[1] - s.jvals[0]) : 1.0;

  auto nearestI = [&](double x) -> int {
    int ii = (s.si > 1) ? static_cast<int>(std::lround((x - s.ivals[0]) / stepI)) : 0;
    return std::max(0, std::min(s.si - 1, ii));
  };
  auto nearestJ = [&](double y) -> int {
    int jj = (s.sj > 1) ? static_cast<int>(std::lround((y - s.jvals[0]) / stepJ)) : 0;
    return std::max(0, std::min(s.sj - 1, jj));
  };

  auto strongAtEdge = [&](const FitPeak2DLocal& from,
                          const FitPeak2DLocal& toward) -> bool {
    const double vx = toward.x - from.x;
    const double vy = toward.y - from.y;
    const double norm = std::sqrt((vx * vx) / (rx * rx) +
                                  (vy * vy) / (ry * ry));
    if (norm <= 1.0) return false;
    const double ex = from.x + vx / norm;
    const double ey = from.y + vy / norm;
    const int ei = nearestI(ex);
    const int ej = nearestJ(ey);
    const int ti = std::max(0, std::min(s.si - 1, toward.i));
    const int tj = std::max(0, std::min(s.sj - 1, toward.j));
    const double edgeI = std::fabs(static_cast<double>(s.DI[ei + ej * s.si]));
    const double peakI = std::fabs(static_cast<double>(s.DI[ti + tj * s.si]));
    return peakI > 0.0 && edgeI > 0.5 * peakI;
  };

  return strongAtEdge(a, b) || strongAtEdge(b, a);
}


// FIT-only area-normalised GLORE evaluator.  Keep spinUnidec's mixing
// convention: g=0 is Gaussian, g=1 is Lorentzian.  The generic Peak()
// routines remain height-normalised and are intentionally untouched.
static inline Peak1DParts EvaluatePeak1DNormalizedFit(double x, double x0,
                                                       double sig, double lor,
                                                       double voigt)
{
  const double ln2 = std::log(2.0);
  const double pi = std::acos(-1.0);
  const double ws = std::max(std::fabs(sig), 1e-12);
  const double wl = std::max(std::fabs(lor), 1e-12);
  const double dx = x - x0;
  const double dx2 = dx * dx;

  const double gauss = (2.0 * std::sqrt(ln2 / pi) / ws) *
                       std::exp(-4.0 * ln2 * dx2 / (ws * ws));
  const double gd_x0 = gauss * (8.0 * ln2 * dx / (ws * ws));
  const double gd_w = gauss * (-1.0 / ws +
                       8.0 * ln2 * dx2 / (ws * ws * ws));

  const double lden = wl * wl + 4.0 * dx2;
  const double lorentz = 2.0 * wl / (pi * lden);
  const double ld_x0 = lorentz * (8.0 * dx / lden);
  const double ld_w = lorentz * (1.0 / wl - 2.0 * wl / lden);

  Peak1DParts out;
  out.value = (1.0 - voigt) * gauss + voigt * lorentz;
  out.d_x0 = (1.0 - voigt) * gd_x0 + voigt * ld_x0;
  out.d_sig = (1.0 - voigt) * gd_w;
  out.d_lor = voigt * ld_w;
  out.d_voigt = lorentz - gauss;
  return out;
}

static inline double BasisValue2D(const FitPeak2DLocal& pk, double x, double y)
{
  return EvaluatePeak1DNormalizedFit(x, pk.x, pk.sig1, pk.lor1, pk.voigt1).value *
         EvaluatePeak1DNormalizedFit(y, pk.y, pk.sig2, pk.lor2, pk.voigt2).value;
}


static void BuildEllipsePixels2D(const slice2D& s,
                               const FitGroup2DLocal& g,
                               const std::vector<FitPeak2DLocal>& peaks,
                               double radIppm,
                               double radJppm,
                               std::vector<int>& pixels)
{
  pixels.clear();
  pixels.reserve((g.maxI - g.minI + 1) * (g.maxJ - g.minJ + 1));

  for (int j = g.minJ; j <= g.maxJ; ++j)
  {
    const double y = s.jvals[j];
    for (int i = g.minI; i <= g.maxI; ++i)
    {
      const double x = s.ivals[i];

      bool keep = false;
      for (int pm : g.members)
      {
        const FitPeak2DLocal& pk = peaks[pm];
        if (InsideEllipse2D(x, y, pk.x, pk.y, radIppm, radJppm))
        {
          keep = true;
          break;
        }
      }

      if (keep)
        pixels.push_back(i + j * s.si);
    }
  }
}

static double FitGroupAtFixedShape2D(const slice2D& s,
                                   const FitGroup2DLocal& g,
                                   const std::vector<FitPeak2DLocal>& peaks,
                                   const std::vector<int>& pixels,
                                   std::vector<double>& coeff)
{
  const int npeak = static_cast<int>(g.members.size());
  if (npeak == 0)
  {
    coeff.clear();
    return 0.0;
  }

  std::vector<double> ata(npeak * npeak, 0.0);
  std::vector<double> aty(npeak, 0.0);

  for (int idx : pixels)
  {
    const int i = idx % s.si;
    const int j = idx / s.si;
    const double yy = static_cast<double>(s.DI[idx]);

    std::vector<double> row(npeak, 0.0);
    for (int p = 0; p < npeak; ++p)
      row[p] = BasisValue2D(peaks[g.members[p]], s.ivals[i], s.jvals[j]);

    for (int a = 0; a < npeak; ++a)
    {
      aty[a] += row[a] * yy;
      for (int b = 0; b <= a; ++b)
        ata[a * npeak + b] += row[a] * row[b];
    }
  }

  for (int a = 0; a < npeak; ++a)
    for (int b = 0; b < a; ++b)
      ata[b * npeak + a] = ata[a * npeak + b];

  if (!SolveLinearSystem(ata, aty, coeff, npeak))
    coeff.assign(npeak, 0.0);

  double sse = 0.0;
  for (int idx : pixels)
  {
    const int i = idx % s.si;
    const int j = idx / s.si;
    double pred = 0.0;
    for (int p = 0; p < npeak; ++p)
      pred += coeff[p] * BasisValue2D(peaks[g.members[p]], s.ivals[i], s.jvals[j]);

    const double r = static_cast<double>(s.DI[idx]) - pred;
    sse += r * r;
  }

  return sse;
}

static void PredictGroupToBuffer2D(const slice2D& s,
                                 const FitGroup2DLocal& g,
                                 const std::vector<FitPeak2DLocal>& peaks,
                                 const std::vector<int>& pixels,
                                 const std::vector<double>& coeff,
                                 std::vector<double>& buffer)
{
  const int npeak = static_cast<int>(g.members.size());
  for (int idx : pixels)
  {
    const int i = idx % s.si;
    const int j = idx / s.si;
    double pred = 0.0;
    for (int p = 0; p < npeak; ++p)
      pred += coeff[p] * BasisValue2D(peaks[g.members[p]], s.ivals[i], s.jvals[j]);
    buffer[idx] = pred;
  }
}

void slice2D::BuildPeakListFromDB(double threshold,
                                  std::vector<FitPeak2DLocal>& peaks)
{
  peaks.clear();

  // Match the numeric names written by slice2D::correlate() for an
  // unrestricted 2D deconvolution (i-major, then j-major traversal).
  std::vector<int> peakNumber(size, 0);
  int ordinal = 0;
  for (int i = 0; i < si; ++i)
    for (int j = 0; j < sj; ++j)
    {
      const int ii = i + j * si;
      if (std::fabs(static_cast<double>(DB[ii])) > threshold)
        peakNumber[ii] = ++ordinal;
    }

  for (int j = 0; j < sj; ++j)
  {
    for (int i = 0; i < si; ++i)
    {
      const int ii = i + j * si;
      const double v = static_cast<double>(DB[ii]);
      if (std::fabs(v) <= threshold)
        continue;

      FitPeak2DLocal p;
      p.i = i;
      p.j = j;
      p.ii = ii;
      p.x = ivals[i];
      p.y = jvals[j];
      p.raw = v;
      p.intensity = v;
      p.fitted = 0.0;
      p.name = std::to_string(peakNumber[ii]);
      p.sig1 = sig1;
      p.sig2 = sig2;
      p.lor1 = lor1;
      p.lor2 = lor2;
      p.voigt1 = voigt1;
      p.voigt2 = voigt2;
      p.group = -1;
      peaks.push_back(p);
    }
  }
}

void slice2D::BuildPeakListFromReference(std::vector<FitPeak2DLocal>& peaks)
{
  peaks.clear();
  for (size_t n = 0; n < peakList.size(); ++n)
  {
    const peakEntry& ref = peakList[n];
    const int i = DoIndex(ref.x, ivals, si);
    const int j = DoIndex(ref.y, jvals, sj);
    if (i < 0 || i >= si || j < 0 || j >= sj)
      continue;

    const int ii = i + j * si;
    FitPeak2DLocal p;
    p.i = i;
    p.j = j;
    p.ii = ii;
    p.x = ref.x;
    p.y = ref.y;
    p.raw = static_cast<double>(DB[ii]);
    p.intensity = p.raw;
    if (std::fabs(p.intensity) < 1.0e-30)
      p.intensity = static_cast<double>(DI[ii]);
    p.fitted = 0.0;
    p.name = ref.name;
    p.sig1 = sig1;
    p.sig2 = sig2;
    p.lor1 = lor1;
    p.lor2 = lor2;
    p.voigt1 = voigt1;
    p.voigt2 = voigt2;
    p.group = -1;
    peaks.push_back(p);
  }
}

std::vector<FitGroup2DLocal> slice2D::BuildGroups(
                                                const std::vector<FitPeak2DLocal>& peaks,
                                                double radIppm,
                                                double radJppm)
{
  std::vector<FitGroup2DLocal> groups;
  const int n = static_cast<int>(peaks.size());
  if (n == 0)
    return groups;

  const double di = (si > 1) ? std::fabs(ivals[1] - ivals[0]) : 1.0;
  const double dj = (sj > 1) ? std::fabs(jvals[1] - jvals[0]) : 1.0;
  const double rx = std::max(radIppm, di);
  const double ry = std::max(radJppm, dj);
  const int cellX = std::max(1, static_cast<int>(std::floor(rx / di)));
  const int cellY = std::max(1, static_cast<int>(std::floor(ry / dj)));

  auto keyOf = [](int a, int b) -> long long
  {
    return (static_cast<long long>(a) << 32) ^ static_cast<unsigned int>(b);
  };

  std::unordered_map<long long, std::vector<int> > bins;
  bins.reserve(n * 2);

  for (int p = 0; p < n; ++p)
    bins[keyOf(peaks[p].i / cellX, peaks[p].j / cellY)].push_back(p);

  std::vector<int> parent(n);
  std::iota(parent.begin(), parent.end(), 0);

  auto findp = [&](int x) -> int
  {
    while (parent[x] != x)
    {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };

  auto unite = [&](int a, int b)
  {
    a = findp(a);
    b = findp(b);
    if (a != b)
      parent[b] = a;
  };

  for (int p = 0; p < n; ++p)
  {
    const int bx = peaks[p].i / cellX;
    const int by = peaks[p].j / cellY;

    for (int dx = -4; dx <= 4; ++dx)
    {
      for (int dy = -4; dy <= 4; ++dy)
      {
        const auto it = bins.find(keyOf(bx + dx, by + dy));
        if (it == bins.end())
          continue;

        for (int q : it->second)
        {
          if (p == q) continue;
          const bool geometricOverlap =
              EllipseOverlap2D(peaks[p].x, peaks[p].y,
                               peaks[q].x, peaks[q].y,
                               radIppm, radJppm);
          const bool intensityOverlap = !geometricOverlap &&
              EdgeIntensityOverlap2D(*this, peaks[p], peaks[q],
                                     radIppm, radJppm);
          if (geometricOverlap || intensityOverlap)
          {
            if (intensityOverlap && findp(p) != findp(q))
              std::cout << "FIT grouping: merging peaks " << peaks[p].name
                        << " and " << peaks[q].name
                        << " because DI remains high at the ROI edge" << std::endl;
            unite(p, q);
          }
        }
      }
    }
  }

  std::unordered_map<int, int> rootToGroup;
  for (int p = 0; p < n; ++p)
  {
    const int root = findp(p);
    int gid;
    auto it = rootToGroup.find(root);
    if (it == rootToGroup.end())
    {
      gid = static_cast<int>(groups.size());
      rootToGroup[root] = gid;
      groups.push_back(FitGroup2DLocal());
    }
    else
    {
      gid = it->second;
    }
    groups[gid].members.push_back(p);
  }

  const int padI = std::max(1, static_cast<int>(std::ceil((radIppm + 4.0 * std::max(sig1, lor1)) / di)));
  const int padJ = std::max(1, static_cast<int>(std::ceil((radJppm + 4.0 * std::max(sig2, lor2)) / dj)));

  for (auto& g : groups)
  {
    int minI = si - 1;
    int maxI = 0;
    int minJ = sj - 1;
    int maxJ = 0;

    for (int p : g.members)
    {
      minI = std::min(minI, peaks[p].i);
      maxI = std::max(maxI, peaks[p].i);
      minJ = std::min(minJ, peaks[p].j);
      maxJ = std::max(maxJ, peaks[p].j);
    }

    g.minI = std::max(0, minI - padI);
    g.maxI = std::min(si - 1, maxI + padI);
    g.minJ = std::max(0, minJ - padJ);
    g.maxJ = std::min(sj - 1, maxJ + padJ);
  }

  return groups;
}

// Reconstruct the authoritative FIT groups from the group IDs assigned when
// the common 2D lineshape fit was created.  Do not call BuildGroups again
// after fitting: fitted peak centres can move enough to change geometric
// overlap and hence renumber/repartition groups.  Protocol3P must keep the
// exact same membership for slice intensities and FUDA-style .dat output.
static bool BuildStoredFitGroups2D(slice2D& s,
                                   const std::vector<FitPeak2DLocal>& peaks,
                                   double radIppm, double radJppm,
                                   std::vector<FitGroup2DLocal>& groups)
{
  groups.clear();
  if (peaks.empty()) return true;

  int maxGroup = -1;
  for (size_t p = 0; p < peaks.size(); ++p)
  {
    if (peaks[p].group < 0)
    {
      std::cerr << "FIT output: peak " << peaks[p].name
                << " has no stored overlap group" << std::endl;
      return false;
    }
    maxGroup = std::max(maxGroup, peaks[p].group);
  }
  groups.resize(maxGroup + 1);
  for (size_t p = 0; p < peaks.size(); ++p)
    groups[peaks[p].group].members.push_back(static_cast<int>(p));

  const double di = (s.si > 1) ? std::fabs(s.ivals[1] - s.ivals[0]) : 1.0;
  const double dj = (s.sj > 1) ? std::fabs(s.jvals[1] - s.jvals[0]) : 1.0;
  const int padI = std::max(1, static_cast<int>(std::ceil(radIppm / di)) + 1);
  const int padJ = std::max(1, static_cast<int>(std::ceil(radJppm / dj)) + 1);

  for (size_t g = 0; g < groups.size(); ++g)
  {
    if (groups[g].members.empty())
    {
      std::cerr << "FIT output: stored overlap group " << g
                << " has no members" << std::endl;
      return false;
    }
    int minI = s.si - 1, maxI = 0, minJ = s.sj - 1, maxJ = 0;
    for (int pm : groups[g].members)
    {
      // Bounds follow the final fitted centre, but membership remains fixed.
      const int i = s.DoIndex(peaks[pm].x, s.ivals, s.si);
      const int j = s.DoIndex(peaks[pm].y, s.jvals, s.sj);
      minI = std::min(minI, i); maxI = std::max(maxI, i);
      minJ = std::min(minJ, j); maxJ = std::max(maxJ, j);
    }
    groups[g].minI = std::max(0, minI - padI);
    groups[g].maxI = std::min(s.si - 1, maxI + padI);
    groups[g].minJ = std::max(0, minJ - padJ);
    groups[g].maxJ = std::min(s.sj - 1, maxJ + padJ);
  }
  return true;
}

static constexpr int kFitOptimisationGuessAndCheck = 0;
static constexpr int kFitOptimisationLevenbergMarquardt = 1;
static constexpr int kFitOptimisationGuessThenLM = 2;


struct PeakShapeParams2D
{
  double intensity;
  double x;
  double y;
  double sig1;
  double sig2;
  double lor1;
  double lor2;
  double voigt1;
  double voigt2;
};

static inline PeakShapeParams2D MakeParamsFromPeak(const FitPeak2DLocal& pk)
{
  PeakShapeParams2D p;
  p.intensity = pk.intensity;
  p.x = pk.x;
  p.y = pk.y;
  p.sig1 = pk.sig1;
  p.sig2 = pk.sig2;
  p.lor1 = pk.lor1;
  p.lor2 = pk.lor2;
  p.voigt1 = pk.voigt1;
  p.voigt2 = pk.voigt2;
  return p;
}

static inline void ApplyParamsToPeak2D(FitPeak2DLocal& pk, const PeakShapeParams2D& p)
{
  pk.intensity = p.intensity;
  pk.fitted = p.intensity;
  pk.x = p.x;
  pk.y = p.y;
  pk.sig1 = p.sig1;
  pk.sig2 = p.sig2;
  pk.lor1 = p.lor1;
  pk.lor2 = p.lor2;
  pk.voigt1 = p.voigt1;
  pk.voigt2 = p.voigt2;
}


static inline double PeakValueAndDerivatives2D(const PeakShapeParams2D& p,
                                             double x, double y,
                                             std::array<double, 9>& deriv)
{
  const Peak1DParts px = EvaluatePeak1DNormalizedFit(x, p.x, p.sig1, p.lor1, p.voigt1);
  const Peak1DParts py = EvaluatePeak1DNormalizedFit(y, p.y, p.sig2, p.lor2, p.voigt2);

  const double basis = px.value * py.value;
  const double model = p.intensity * basis;

  deriv[0] = basis;
  deriv[1] = p.intensity * px.d_x0 * py.value;
  deriv[2] = p.intensity * px.d_sig * py.value;
  deriv[3] = p.intensity * px.d_lor * py.value;
  deriv[4] = p.intensity * px.d_voigt * py.value;
  deriv[5] = p.intensity * px.value * py.d_x0;
  deriv[6] = p.intensity * px.value * py.d_sig;
  deriv[7] = p.intensity * px.value * py.d_lor;
  deriv[8] = p.intensity * px.value * py.d_voigt;

  return model;
}

static inline void ClampParams2D(PeakShapeParams2D& p,
                               const slice2D& s)
{
  const double xmin = std::min(s.ivals[0], s.ivals[s.si - 1]);
  const double xmax = std::max(s.ivals[0], s.ivals[s.si - 1]);
  const double ymin = std::min(s.jvals[0], s.jvals[s.sj - 1]);
  const double ymax = std::max(s.jvals[0], s.jvals[s.sj - 1]);

  if (!std::isfinite(p.intensity))
    p.intensity = 0.0;
  if (!std::isfinite(p.x))
    p.x = 0.5 * (xmin + xmax);
  if (!std::isfinite(p.y))
    p.y = 0.5 * (ymin + ymax);
  if (!std::isfinite(p.sig1) || p.sig1 <= 1e-12)
    p.sig1 = 1e-12;
  if (!std::isfinite(p.sig2) || p.sig2 <= 1e-12)
    p.sig2 = 1e-12;
  if (!std::isfinite(p.lor1) || p.lor1 <= 1e-12)
    p.lor1 = 1e-12;
  if (!std::isfinite(p.lor2) || p.lor2 <= 1e-12)
    p.lor2 = 1e-12;
  if (!std::isfinite(p.voigt1))
    p.voigt1 = 0.0;
  if (!std::isfinite(p.voigt2))
    p.voigt2 = 0.0;

  p.x = clampd(p.x, xmin, xmax);
  p.y = clampd(p.y, ymin, ymax);
  p.voigt1 = clampd(p.voigt1, 0.0, 1.0);
  p.voigt2 = clampd(p.voigt2, 0.0, 1.0);
}

// LM refinement is deliberately local: the deconvolution/guess-and-check
// stage has already located the peak.  Permit only a half-ROI movement
// around that starting position while retaining the existing physical
// linewidth and Voigt constraints.
static inline void ClampParams2DLocal(PeakShapeParams2D& p,
                                      const slice2D& s,
                                      double x0, double y0,
                                      double halfRadIppm,
                                      double halfRadJppm,
                                      bool restrictWidths = false,
                                      double sig1Ref = 0.0, double sig2Ref = 0.0,
                                      double lor1Ref = 0.0, double lor2Ref = 0.0)
{
  ClampParams2D(p, s);

  const double xminSpec = std::min(s.ivals[0], s.ivals[s.si - 1]);
  const double xmaxSpec = std::max(s.ivals[0], s.ivals[s.si - 1]);
  const double yminSpec = std::min(s.jvals[0], s.jvals[s.sj - 1]);
  const double ymaxSpec = std::max(s.jvals[0], s.jvals[s.sj - 1]);

  const double xmin = std::max(xminSpec, x0 - std::fabs(halfRadIppm));
  const double xmax = std::min(xmaxSpec, x0 + std::fabs(halfRadIppm));
  const double ymin = std::max(yminSpec, y0 - std::fabs(halfRadJppm));
  const double ymax = std::min(ymaxSpec, y0 + std::fabs(halfRadJppm));

  p.x = clampd(p.x, xmin, xmax);
  p.y = clampd(p.y, ymin, ymax);

  if (restrictWidths)
  {
    // FUDA-like mode: one exact FWHM per dimension.  The Gaussian and
    // Lorentzian components share that width; g alone controls the mixture.
    const double w1 = std::max(1e-12, p.sig1);
    const double w2 = std::max(1e-12, p.sig2);
    p.sig1 = p.lor1 = w1;
    p.sig2 = p.lor2 = w2;
  }
}

static inline int ParamOffset(int peakIndex)
{
  return peakIndex * 9;
}

static double EvaluateGroupLM(const slice2D& s,
                              const FitGroup2DLocal& g,
                              const std::vector<int>& pixels,
                              const std::vector<PeakShapeParams2D>& params,
                              std::vector<double>* model,
                              std::vector<double>* raw,
                              std::vector<double>* jtj,
                              std::vector<double>* jtr)
{
  const int npeak = static_cast<int>(params.size());
  const int nparam = npeak * 9;

  double sse = 0.0;
  for (int idx : pixels)
  {
    const int i = idx % s.si;
    const int j = idx / s.si;
    const double yy = static_cast<double>(s.DI[idx]);

    std::vector<std::array<double, 9> > derivs(npeak);
    double pred = 0.0;
    for (int p = 0; p < npeak; ++p)
      pred += PeakValueAndDerivatives2D(params[p], s.ivals[i], s.jvals[j], derivs[p]);

    const double residual = yy - pred;
    sse += residual * residual;

    if (model)
      (*model)[idx] = pred;
    if (raw)
      (*raw)[idx] = yy;

    if (jtj && jtr)
    {
      // Accumulate the complete normal equations.  This must include
      // cross-peak terms as well as the 9x9 block for each individual peak;
      // otherwise an overlapped group is not a true simultaneous LM fit.
      for (int p = 0; p < npeak; ++p)
      {
        const int offP = ParamOffset(p);
        const std::array<double, 9>& dp = derivs[p];
        for (int a = 0; a < 9; ++a)
        {
          const int ia = offP + a;
          (*jtr)[ia] += dp[a] * residual;
          for (int q = 0; q <= p; ++q)
          {
            const int offQ = ParamOffset(q);
            const std::array<double, 9>& dq = derivs[q];
            const int bmax = (q == p) ? a : 8;
            for (int b = 0; b <= bmax; ++b)
            {
              const int ib = offQ + b;
              (*jtj)[ia * nparam + ib] += dp[a] * dq[b];
            }
          }
        }
      }
    }
  }

  if (jtj && jtr)
  {
    for (int a = 0; a < nparam; ++a)
      for (int b = 0; b < a; ++b)
        (*jtj)[b * nparam + a] = (*jtj)[a * nparam + b];
  }

  return sse;
}


static double FitOneGroup2DGuess(const slice2D& s,
                                 const FitGroup2DLocal& g,
                                 std::vector<FitPeak2DLocal>& peaks,
                                 int maxIter,
                                 std::vector<double>& model,
                                 std::vector<double>& raw,
                                 double radIppm,
                                 double radJppm,
                                 bool restrictWidths)
{
  const int npeak = static_cast<int>(g.members.size());
  if (npeak == 0)
    return 0.0;

  // In FUDA-like shared-width mode the Gaussian FWHM is the single
  // starting linewidth for each dimension and the Lorentzian FWHM is tied
  // to it exactly from the first guess/check evaluation onward.
  if (restrictWidths)
  {
    for (int member : g.members)
    {
      peaks[member].lor1 = peaks[member].sig1;
      peaks[member].lor2 = peaks[member].sig2;
    }
  }

  std::vector<int> pixels;
  BuildEllipsePixels2D(s, g, peaks, radIppm, radJppm, pixels);

  std::vector<double> coeff(npeak, 0.0);
  double best = FitGroupAtFixedShape2D(s, g, peaks, pixels, coeff);

  auto probeOneParam = [&](int memberIndex, int paramIndex)
  {
    FitPeak2DLocal& pk = peaks[g.members[memberIndex]];
    double* target = 0;
    bool sharedF1 = false, sharedF2 = false;
    switch (paramIndex)
    {
      case 0: target = &pk.sig1; sharedF1 = restrictWidths; break;
      case 1: target = &pk.sig2; sharedF2 = restrictWidths; break;
      case 2: if (restrictWidths) return; target = &pk.lor1; break;
      case 3: if (restrictWidths) return; target = &pk.lor2; break;
      case 4: target = &pk.voigt1; break;
      case 5: target = &pk.voigt2; break;
      default: return;
    }

    const double old = *target;
    double candidates[5];
    int nc = 0;

    if (paramIndex <= 3)
    {
      candidates[nc++] = std::max(1e-12, old * 0.75);
      candidates[nc++] = std::max(1e-12, old * 0.875);
      candidates[nc++] = old;
      candidates[nc++] = old * 1.125;
      candidates[nc++] = old * 1.25;
    }
    else
    {
      candidates[nc++] = clampd(old - 0.20, 0.0, 1.0);
      candidates[nc++] = clampd(old - 0.10, 0.0, 1.0);
      candidates[nc++] = old;
      candidates[nc++] = clampd(old + 0.10, 0.0, 1.0);
      candidates[nc++] = clampd(old + 0.20, 0.0, 1.0);
    }

    double localBest = best;
    double localParam = old;
    std::vector<double> trialCoeff = coeff;
    std::vector<double> bestCoeffLocal = coeff;

    for (int c = 0; c < nc; ++c)
    {
      *target = candidates[c];
      if (sharedF1) pk.lor1 = candidates[c];
      if (sharedF2) pk.lor2 = candidates[c];
      trialCoeff = coeff;
      const double sse = FitGroupAtFixedShape2D(s, g, peaks, pixels, trialCoeff);
      if (sse < localBest)
      {
        localBest = sse;
        localParam = candidates[c];
        bestCoeffLocal = trialCoeff;
      }
    }

    *target = localParam;
    if (sharedF1) pk.lor1 = localParam;
    if (sharedF2) pk.lor2 = localParam;
    coeff = bestCoeffLocal;
    best = localBest;
  };

  for (int outer = 0; outer < maxIter; ++outer)
  {
    const double prev = best;

    for (int p = 0; p < npeak; ++p)
    {
      for (int m = 0; m < 6; ++m)
        probeOneParam(p, m);
    }

    if (std::fabs(prev - best) / (std::fabs(prev) + 1e-12) < 1e-6)
      break;
  }

  best = FitGroupAtFixedShape2D(s, g, peaks, pixels, coeff);

  for (int p = 0; p < npeak; ++p)
  {
    FitPeak2DLocal& pk = peaks[g.members[p]];
    pk.intensity = coeff[p];
    pk.fitted = pk.intensity;
  }

  PredictGroupToBuffer2D(s, g, peaks, pixels, coeff, model);

  for (int idx : pixels)
    raw[idx] = static_cast<double>(s.DI[idx]);

  return best;
}


// FUDA-like shared-width LM model used when FitWidthRestrict is enabled.
// Each dimension has one FWHM shared exactly by its Gaussian and Lorentzian
// components: [I, x, w1, g1, y, w2, g2] per peak.
static double EvaluateGroupLMSharedWidth(const slice2D& s,
                              const FitGroup2DLocal& g,
                              const std::vector<int>& pixels,
                              const std::vector<PeakShapeParams2D>& params,
                              std::vector<double>* model,
                              std::vector<double>* raw,
                              std::vector<double>* jtj,
                              std::vector<double>* jtr)
{
  const int npeak = static_cast<int>(params.size());
  const int nper = 7;
  const int nparam = npeak * nper;
  double sse = 0.0;
  for (int idx : pixels)
  {
    const int i = idx % s.si, j = idx / s.si;
    const double yy = static_cast<double>(s.DI[idx]);
    std::vector<std::array<double,7> > derivs(npeak);
    double pred = 0.0;
    for (int p = 0; p < npeak; ++p)
    {
      const PeakShapeParams2D& pp = params[p];
      const Peak1DParts px = EvaluatePeak1DNormalizedFit(s.ivals[i], pp.x, pp.sig1, pp.sig1, pp.voigt1);
      const Peak1DParts py = EvaluatePeak1DNormalizedFit(s.jvals[j], pp.y, pp.sig2, pp.sig2, pp.voigt2);
      const double basis = px.value * py.value;
      pred += pp.intensity * basis;
      derivs[p][0] = basis;
      derivs[p][1] = pp.intensity * px.d_x0 * py.value;
      // Because w is shared, d/dw = d/dsig + d/dlor.
      derivs[p][2] = pp.intensity * (px.d_sig + px.d_lor) * py.value;
      derivs[p][3] = pp.intensity * px.d_voigt * py.value;
      derivs[p][4] = pp.intensity * px.value * py.d_x0;
      derivs[p][5] = pp.intensity * px.value * (py.d_sig + py.d_lor);
      derivs[p][6] = pp.intensity * px.value * py.d_voigt;
    }
    const double residual = yy - pred;
    sse += residual * residual;
    if (model) (*model)[idx] = pred;
    if (raw) (*raw)[idx] = yy;
    if (jtj && jtr)
    {
      for (int p=0;p<npeak;++p) for(int a=0;a<nper;++a)
      {
        const int ia=p*nper+a; (*jtr)[ia]+=derivs[p][a]*residual;
        for(int q=0;q<=p;++q)
        {
          const int bmax=(q==p)?a:nper-1;
          for(int b=0;b<=bmax;++b)
            (*jtj)[ia*nparam + q*nper+b] += derivs[p][a]*derivs[q][b];
        }
      }
    }
  }
  // Penalise shared linewidths below the digital sampling interval.  The
  // residual scale is the RMS signal in this fit region, so the penalty is
  // expressed in the same units as the spectral residuals and adapts to data.
  const double minW1 = (s.si > 1) ? std::fabs(s.ivals[1] - s.ivals[0]) : 1e-12;
  const double minW2 = (s.sj > 1) ? std::fabs(s.jvals[1] - s.jvals[0]) : 1e-12;
  double rms = 0.0;
  for (int idx : pixels) { const double v = static_cast<double>(s.DI[idx]); rms += v*v; }
  rms = pixels.empty() ? 1.0 : std::sqrt(rms / static_cast<double>(pixels.size()));
  if (!(rms > 0.0) || !std::isfinite(rms)) rms = 1.0;
  for (int p=0; p<npeak; ++p)
  {
    const double widths[2] = {params[p].sig1, params[p].sig2};
    const double mins[2] = {minW1, minW2};
    const int slots[2] = {2, 5};
    for (int d=0; d<2; ++d)
    {
      if (widths[d] >= mins[d]) continue;
      const double rpen = rms * (mins[d] - widths[d]) / mins[d];
      const double jpen = rms / mins[d];
      sse += rpen * rpen;
      if (jtj && jtr)
      {
        const int a = p*nper + slots[d];
        (*jtr)[a] += jpen * rpen;
        (*jtj)[a*nparam+a] += jpen * jpen;
      }
    }
  }

  if(jtj&&jtr) for(int a=0;a<nparam;++a) for(int b=0;b<a;++b)
    (*jtj)[b*nparam+a]=(*jtj)[a*nparam+b];
  return sse;
}

static double FitOneGroup2DLM(const slice2D& s,
                              const FitGroup2DLocal& g,
                              std::vector<FitPeak2DLocal>& peaks,
                              int maxIter,
                              std::vector<double>& model,
                              std::vector<double>& raw,
                              double radIppm,
                              double radJppm,
                              bool restrictWidths)
{
  const int npeak = static_cast<int>(g.members.size());
  if (npeak == 0)
    return 0.0;

  std::vector<int> pixels;
  BuildEllipsePixels2D(s, g, peaks, radIppm, radJppm, pixels);
  if (pixels.empty())
    return 0.0;

  std::vector<double> coeff(npeak, 0.0);
  const double startSSE = FitGroupAtFixedShape2D(s, g, peaks, pixels, coeff);

  // The input peaks are the LM starting solution.  Record their positions
  // before refinement so each centre can move by at most half an extraction
  // radius in either dimension.
  std::vector<double> startX(npeak, 0.0);
  std::vector<double> startY(npeak, 0.0);
  std::vector<double> startSig1(npeak, 0.0), startSig2(npeak, 0.0);
  std::vector<double> startLor1(npeak, 0.0), startLor2(npeak, 0.0);
  for (int p = 0; p < npeak; ++p)
  {
    const FitPeak2DLocal& startPk = peaks[g.members[p]];
    startX[p] = startPk.x;
    startY[p] = startPk.y;
    startSig1[p] = startPk.sig1;
    startSig2[p] = startPk.sig2;
    startLor1[p] = startPk.lor1;
    startLor2[p] = startPk.lor2;
  }

  std::vector<PeakShapeParams2D> params(npeak);
  for (int p = 0; p < npeak; ++p)
  {
    const FitPeak2DLocal& pk = peaks[g.members[p]];
    params[p] = MakeParamsFromPeak(pk);
    if (restrictWidths)
    {
      const double w1 = 0.5 * (params[p].sig1 + params[p].lor1);
      const double w2 = 0.5 * (params[p].sig2 + params[p].lor2);
      params[p].sig1 = params[p].lor1 = w1;
      params[p].sig2 = params[p].lor2 = w2;
    }
    params[p].intensity = coeff[p];
    if (!std::isfinite(params[p].intensity))
      params[p].intensity = pk.raw;
    ClampParams2DLocal(params[p], s, startX[p], startY[p],
                       0.5 * radIppm, 0.5 * radJppm, restrictWidths,
                       startSig1[p], startSig2[p], startLor1[p], startLor2[p]);
  }

  std::vector<PeakShapeParams2D> bestParams = params;
  double bestSSE = startSSE;
  double prevAcceptedSSE = startSSE;

  const int nper = restrictWidths ? 7 : 9;
  const int nparam = npeak * nper;
  const int maxAttempts = 8;
  double lambda = 1e-3;

  for (int iter = 0; iter < maxIter; ++iter)
  {
    bool accepted = false;

    for (int attempt = 0; attempt < maxAttempts; ++attempt)
    {
      std::vector<double> jtj(nparam * nparam, 0.0);
      std::vector<double> jtr(nparam, 0.0);

      const double currentSSE = restrictWidths ?
        EvaluateGroupLMSharedWidth(s, g, pixels, params, 0, 0, &jtj, &jtr) :
        EvaluateGroupLM(s, g, pixels, params, 0, 0, &jtj, &jtr);

      if (!std::isfinite(currentSSE))
      {
        lambda *= 10.0;
        continue;
      }

      for (int p = 0; p < nparam; ++p)
      {
        const double diag = jtj[p * nparam + p];
        jtj[p * nparam + p] += lambda * std::max(1.0, diag);
      }

      std::vector<double> delta;
      if (!SolveLinearSystem(jtj, jtr, delta, nparam))
      {
        lambda *= 10.0;
        continue;
      }

      std::vector<PeakShapeParams2D> trial = params;
      for (int p = 0; p < npeak; ++p)
      {
        if (restrictWidths)
        {
          const int off = p * 7;
          trial[p].intensity += delta[off + 0];
          trial[p].x         += delta[off + 1];
          trial[p].sig1      += delta[off + 2];
          trial[p].lor1       = trial[p].sig1;
          trial[p].voigt1    += delta[off + 3];
          trial[p].y         += delta[off + 4];
          trial[p].sig2      += delta[off + 5];
          trial[p].lor2       = trial[p].sig2;
          trial[p].voigt2    += delta[off + 6];
        }
        else
        {
          const int off = ParamOffset(p);
          trial[p].intensity += delta[off + 0];
          trial[p].x         += delta[off + 1];
          trial[p].sig1      += delta[off + 2];
          trial[p].lor1      += delta[off + 3];
          trial[p].voigt1    += delta[off + 4];
          trial[p].y         += delta[off + 5];
          trial[p].sig2      += delta[off + 6];
          trial[p].lor2      += delta[off + 7];
          trial[p].voigt2    += delta[off + 8];
        }
        ClampParams2DLocal(trial[p], s, startX[p], startY[p],
                           0.5 * radIppm, 0.5 * radJppm, restrictWidths,
                           startSig1[p], startSig2[p], startLor1[p], startLor2[p]);
      }

      const double trialSSE = restrictWidths ?
        EvaluateGroupLMSharedWidth(s, g, pixels, trial, 0, 0, 0, 0) :
        EvaluateGroupLM(s, g, pixels, trial, 0, 0, 0, 0);
      if (std::isfinite(trialSSE) && trialSSE < bestSSE)
      {
        params = trial;
        bestSSE = trialSSE;
        bestParams = trial;
        lambda = std::max(lambda * 0.3, 1e-12);
        accepted = true;
        break;
      }

      lambda *= 10.0;
    }

    if (!accepted)
      break;

    if (std::fabs(prevAcceptedSSE - bestSSE) / (std::fabs(prevAcceptedSSE) + 1e-12) < 1e-10)
      break;

    prevAcceptedSSE = bestSSE;
  }

  for (int p = 0; p < npeak; ++p)
    ApplyParamsToPeak2D(peaks[g.members[p]], bestParams[p]);

  if (restrictWidths)
    EvaluateGroupLMSharedWidth(s, g, pixels, bestParams, &model, &raw, 0, 0);
  else
    EvaluateGroupLM(s, g, pixels, bestParams, &model, &raw, 0, 0);

  return bestSSE;
}

static double FitOneGroup2DGuessThenLM(const slice2D& s,
                                        const FitGroup2DLocal& g,
                                        std::vector<FitPeak2DLocal>& peaks,
                                        int maxIter,
                                        std::vector<double>& model,
                                        std::vector<double>& raw,
                                        double radIppm,
                                        double radJppm,
                                        bool restrictWidths)
{
  // Stage 1: retain the robust historical discrete optimisation.
  const double guessSSE = FitOneGroup2DGuess(s, g, peaks, maxIter,
                                             model, raw, radIppm, radJppm, restrictWidths);

  // Stage 2: refine a copy continuously.  Never allow a failed or poorer LM
  // result to overwrite the established guess-and-check solution.
  std::vector<FitPeak2DLocal> lmPeaks = peaks;
  // Preserve contributions already written by earlier groups; LM only
  // overwrites the pixels belonging to the current group.
  std::vector<double> lmModel = model;
  std::vector<double> lmRaw = raw;
  const double lmSSE = FitOneGroup2DLM(s, g, lmPeaks, maxIter,
                                      lmModel, lmRaw, radIppm, radJppm, restrictWidths);

  const double improveTol = 1e-10;
  const bool acceptLM = std::isfinite(lmSSE) &&
                        (!std::isfinite(guessSSE) ||
                         lmSSE < guessSSE * (1.0 - improveTol));

  std::cout << "FIT group";
  for (size_t m = 0; m < g.members.size(); ++m)
    std::cout << (m == 0 ? " " : ",") << peaks[g.members[m]].name;
  std::cout << ": Guess SSE=" << guessSSE
            << " LM SSE=" << lmSSE
            << (acceptLM ? " LM accepted" : " LM rejected")
            << std::endl;

  if (acceptLM)
  {
    for (int member : g.members)
      peaks[member] = lmPeaks[member];
    model.swap(lmModel);
    raw.swap(lmRaw);
    return lmSSE;
  }

  return guessSSE;
}

static double FitOneGroup2DMode(const slice2D& s,
                                const FitGroup2DLocal& g,
                                std::vector<FitPeak2DLocal>& peaks,
                                int maxIter,
                                std::vector<double>& model,
                                std::vector<double>& raw,
                                double radIppm,
                                double radJppm,
                                int optimisationMode,
                                bool restrictWidths)
{
  if (optimisationMode == kFitOptimisationLevenbergMarquardt)
    return FitOneGroup2DLM(s, g, peaks, maxIter, model, raw, radIppm, radJppm, restrictWidths);
  if (optimisationMode == kFitOptimisationGuessThenLM)
    return FitOneGroup2DGuessThenLM(s, g, peaks, maxIter, model, raw, radIppm, radJppm, restrictWidths);

  return FitOneGroup2DGuess(s, g, peaks, maxIter, model, raw, radIppm, radJppm, restrictWidths);
}

static void FitPeaks2DImpl(slice2D& s,
                           double radIppm,
                           double radJppm,
                           const std::string& paramOut,
                           const std::string& gnuplotOut,
                           int maxIter,
                           double threshold,
                           int optimisationMode,
                           bool useReferencePeakList,
                           bool restrictWidths)
{
  std::vector<FitPeak2DLocal> fitPeaks;
  if (useReferencePeakList)
    s.BuildPeakListFromReference(fitPeaks);
  else
    s.BuildPeakListFromDB(threshold, fitPeaks);

  if (fitPeaks.empty())
  {
    std::cout << "FitPeaks2D: no non-zero DB entries found" << std::endl;
    FILE* fp = fopen(paramOut.c_str(), "w");
    if (fp) fclose(fp);
    fp = fopen(gnuplotOut.c_str(), "w");
    if (fp) fclose(fp);
    return;
  }

  std::vector<FitGroup2DLocal> groups = s.BuildGroups(fitPeaks, radIppm, radJppm);

  for (size_t g = 0; g < groups.size(); ++g)
    for (int idx : groups[g].members)
      fitPeaks[idx].group = static_cast<int>(g);

  std::vector<double> model(s.size, 0.0);
  std::vector<double> raw(s.size, 0.0);
  std::vector<double> groupChi2(groups.size(), 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int g = 0; g < static_cast<int>(groups.size()); ++g)
  {
    std::vector<double> localModel(s.size, 0.0);
    std::vector<double> localRaw(s.size, 0.0);

    groupChi2[g] = FitOneGroup2DMode(s,
                                     groups[g],
                                     fitPeaks,
                                     maxIter,
                                     localModel,
                                     localRaw,
                                     radIppm,
                                     radJppm,
                                     optimisationMode,
                                     restrictWidths);

    for (int j = groups[g].minJ; j <= groups[g].maxJ; ++j)
      for (int i = groups[g].minI; i <= groups[g].maxI; ++i)
      {
        const int ii = i + j * s.si;
        if (localModel[ii] != 0.0)
          model[ii] = localModel[ii];
        if (localRaw[ii] != 0.0)
          raw[ii] = localRaw[ii];
      }
  }

  FILE* pfp = fopen(paramOut.c_str(), "w");
  if (pfp)
  {
    fprintf(pfp, "#peak\txppm\typpm\trawDB\tfittedI\tw1\tw2\tvoigt1\tvoigt2\tgroup\n");
    for (size_t i = 0; i < fitPeaks.size(); ++i)
    {
      const FitPeak2DLocal& p = fitPeaks[i];
      fprintf(pfp, "%zu\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%i\n",
              i, p.x, p.y, p.raw, p.intensity, p.sig1, p.sig2, p.voigt1, p.voigt2, p.group);
    }
    fclose(pfp);
  }

  FILE* gfp = fopen(gnuplotOut.c_str(), "w");
  if (gfp)
  {
    for (size_t g = 0; g < groups.size(); ++g)
    {
      const FitGroup2DLocal& grp = groups[g];
      std::vector<int> pixels;
      BuildEllipsePixels2D(s, grp, fitPeaks, radIppm, radJppm, pixels);

      int lastJ = -999999;
      for (size_t n = 0; n < pixels.size(); ++n)
      {
        const int idx = pixels[n];
        const int i = idx % s.si;
        const int j = idx / s.si;

        if (j != lastJ)
        {
          if (lastJ != -999999)
            fprintf(gfp, "\n");
          lastJ = j;
        }

        double pred = 0.0;
        for (int p = 0; p < static_cast<int>(grp.members.size()); ++p)
          pred += fitPeaks[grp.members[p]].intensity *
                  BasisValue2D(fitPeaks[grp.members[p]], s.ivals[i], s.jvals[j]);

        fprintf(gfp, "%.12e\t%.12e\t%.12e\t%.12e\n",
                s.ivals[i], s.jvals[j], static_cast<double>(s.DI[idx]), pred);
      }
      fprintf(gfp, "\n\n");
    }
    fclose(gfp);
  }

  const double chi2 = std::accumulate(groupChi2.begin(), groupChi2.end(), 0.0);
  std::cout << "FitPeaks2D: peaks=" << fitPeaks.size()
            << " groups=" << groups.size()
            << " chi2=" << chi2 << std::endl;

  s.fittedPeaks2D = fitPeaks;
}

double slice2D::FitOneGroup2D(const FitGroup2DLocal& g,
                             std::vector<FitPeak2DLocal>& peaks,
                             int maxIter,
                             std::vector<double>& model,
                             std::vector<double>& raw,
                             double radIppm,
                             double radJppm)
{
  return FitOneGroup2DGuess(*this, g, peaks, maxIter, model, raw, radIppm, radJppm, false);
}

void slice2D::FitPeaks2D(double radIppm,
                          double radJppm,
                          const std::string& paramOut,
                          const std::string& gnuplotOut,
                          int maxIter,
                          double threshold,
                          bool useReferencePeakList,
                          bool restrictLMWidths)
{
  // Experimental hybrid build: first run the robust historical
  // guess-and-check optimiser, then use that solution to seed a local
  // continuous LM refinement.  LM is accepted only when it lowers SSE.
  int optimisationMode = kFitOptimisationGuessThenLM;
  int iterations = maxIter;
  if (iterations < 0)
    iterations = -iterations;

  FitPeaks2DImpl(*this,
                 radIppm,
                 radJppm,
                 paramOut,
                 gnuplotOut,
                 iterations,
                 threshold,
                 optimisationMode,
                 useReferencePeakList,
                 restrictLMWidths);
}

static std::string SafeFitName2D(const std::string& name, size_t fallback)
{
  std::string out = name.empty() ? std::to_string(fallback + 1) : name;
  for (size_t i = 0; i < out.size(); ++i)
    if (out[i] == '/' || out[i] == '\\') out[i] = '_';
  return out;
}

bool slice2D::WriteFudaFitOutputs(const std::string& fitDir,
                                  double radIppm, double radJppm,
                                  double obs1MHz, double obs2MHz)
{
  if (fittedPeaks2D.empty()) return true;

  if (::mkdir(fitDir.c_str(), 0775) != 0 && errno != EEXIST)
  {
    std::cerr << "Cannot create FIT output directory: " << fitDir
              << " (errno " << errno << ")" << std::endl;
    return false;
  }

  std::vector<FitGroup2DLocal> groups;
  if (!BuildStoredFitGroups2D(*this, fittedPeaks2D, radIppm, radJppm, groups))
    return false;
  const double nanv = std::numeric_limits<double>::quiet_NaN();

  for (size_t n = 0; n < fittedPeaks2D.size(); ++n)
  {
    const FitPeak2DLocal& pk = fittedPeaks2D[n];
    const std::string base = SafeFitName2D(pk.name, n);
    const std::string datName = fitDir + "/" + base + ".dat";
    const std::string outName = fitDir + "/" + base + ".out";

    const int groupIndex = pk.group;
    if (groupIndex < 0 || groupIndex >= static_cast<int>(groups.size()))
    {
      std::cerr << "FIT output: invalid stored group " << groupIndex
                << " for peak " << base << std::endl;
      return false;
    }
    const FitGroup2DLocal& grp = groups[groupIndex];

    FILE* dfp = fopen(datName.c_str(), "w");
    if (!dfp) { std::cerr << "Cannot write FIT data file: " << datName << std::endl; return false; }
    fprintf(dfp, "%11s %11s %11s %11s\n", "# F2(ppm)", "F1(ppm)", "Data", "Calc");
    std::vector<int> pixels;
    BuildEllipsePixels2D(*this, grp, fittedPeaks2D, radIppm, radJppm, pixels);
    int lastJ = -999999;
    for (size_t q = 0; q < pixels.size(); ++q)
    {
      const int idx = pixels[q];
      const int i = idx % si;
      const int j = idx / si;
      if (j != lastJ) { if (lastJ != -999999) fprintf(dfp, "\n"); lastJ = j; }
      double pred = 0.0;
      for (size_t m = 0; m < grp.members.size(); ++m)
      {
        const FitPeak2DLocal& gp = fittedPeaks2D[grp.members[m]];
        pred += gp.intensity * BasisValue2D(gp, ivals[i], jvals[j]);
      }
      fprintf(dfp, "%11.4e %11.4e %11.4e %11.4e\n",
              static_cast<double>(jvals[j]), static_cast<double>(ivals[i]),
              static_cast<double>(DI[idx]), pred);
    }
    fclose(dfp);

    FILE* ofp = fopen(outName.c_str(), "w");
    if (!ofp) { std::cerr << "Cannot write FIT parameter file: " << outName << std::endl; return false; }
    fprintf(ofp, "#\n# %15s%15s\n", "Peak Name", "Overlap_group");
    fprintf(ofp, "# %15s%15d\n#\n", base.c_str(), pk.group);
    fprintf(ofp, "# Input frequencies:\n# %10s%10s%10s\n", "Omega1", "Omega2", "Omega3");
    fprintf(ofp, "# %10.3f%10.3f%10.3f\n", pk.x, pk.y, 0.0);
    fprintf(ofp, "#\n# --------- Results of the fit -------------\n#\n");
    fprintf(ofp, "# Parameter        Value           Esd\n");

    // Shared-width FIT model: Gaussian and Lorentzian components use the
    // same FWHM.  g retains the spinUnidec convention (0=Gaussian, 1=Lorentzian).
    const double w1ppm = pk.sig1;
    const double w2ppm = pk.sig2;
    fprintf(ofp, "  %-12s %14.7e %14.7e\n", "f01(ppm)", pk.x, nanv);
    fprintf(ofp, "  %-12s %14.7e %14.7e\n", "w1(Hz)", w1ppm * obs1MHz, nanv);
    fprintf(ofp, "  %-12s %14.7e %14.7e\n", "g1", pk.voigt1, nanv);
    fprintf(ofp, "  %-12s %14.7e %14.7e\n", "f02(ppm)", pk.y, nanv);
    fprintf(ofp, "  %-12s %14.7e %14.7e\n", "w2(Hz)", w2ppm * obs2MHz, nanv);
    fprintf(ofp, "  %-12s %14.7e %14.7e\n", "g2", pk.voigt2, nanv);
    fprintf(ofp, "  %-12s %14.7e %14.7e\n", "Intensity", pk.intensity, nanv);
    fclose(ofp);
  }
  return true;
}

void slice2D::RebuildDSFromFit(double nsig)
{
  memset(DS, 0, size * SIZEMEM);

  if (fittedPeaks2D.empty())
    return;

  for (size_t p = 0; p < fittedPeaks2D.size(); ++p)
  {
    const FitPeak2DLocal& pk = fittedPeaks2D[p];

    const double wx = nsig * std::max(pk.sig1, pk.lor1);
    const double wy = nsig * std::max(pk.sig2, pk.lor2);

    int i0 = 0, i1 = si - 1;
    int j0 = 0, j1 = sj - 1;

    while (i0 < si && std::fabs(ivals[i0] - pk.x) > wx) ++i0;
    while (i1 >= 0 && std::fabs(ivals[i1] - pk.x) > wx) --i1;
    while (j0 < sj && std::fabs(jvals[j0] - pk.y) > wy) ++j0;
    while (j1 >= 0 && std::fabs(jvals[j1] - pk.y) > wy) --j1;

    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i1 >= si) i1 = si - 1;
    if (j1 >= sj) j1 = sj - 1;

    for (int j = j0; j <= j1; ++j)
    {
      const double gy = EvaluatePeak1DNormalizedFit(jvals[j], pk.y, pk.sig2, pk.lor2, pk.voigt2).value;
      if (gy == 0.0) continue;

      for (int i = i0; i <= i1; ++i)
      {
        const int ii = i + j * si;
        const double gx = EvaluatePeak1DNormalizedFit(ivals[i], pk.x, pk.sig1, pk.lor1, pk.voigt1).value;
        if (gx == 0.0) continue;

        DS[ii] += pk.intensity * gx * gy;
      }
    }
  }
}

#endif

// NOTE: this block intentionally sits after the legacy include guard body in
// older trees; slice2Dfit.cpp is compiled once, so the member definition is
// still available to the linker.
bool slice2D::FitPseudo3DIntensities(const std::vector<double>& stack, int nslices,
                                     double radIppm, double radJppm,
                                     std::vector<std::vector<double> >& intensities,
                                     std::vector<std::vector<double> >& intensityEsd,
                                     const std::string& fitDir,
                                     const std::vector<double>& zvals,
                                     const std::string& zlabel,
                                     double obsI, double obsJ)
{
  if (nslices <= 0 || static_cast<int>(stack.size()) != nslices * size || fittedPeaks2D.empty())
    return false;
  intensities.assign(fittedPeaks2D.size(), std::vector<double>(nslices, 0.0));
  intensityEsd.assign(fittedPeaks2D.size(), std::vector<double>(nslices, 0.0));
  std::vector<FitGroup2DLocal> groups;
  if (!BuildStoredFitGroups2D(*this, fittedPeaks2D, radIppm, radJppm, groups))
    return false;

  for (const FitGroup2DLocal& grp : groups)
  {
    std::vector<int> pixels;
    BuildEllipsePixels2D(*this, grp, fittedPeaks2D, radIppm, radJppm, pixels);
    const int np = static_cast<int>(grp.members.size());
    if (np == 0 || pixels.empty()) continue;

    std::vector<double> ata(np*np, 0.0);
    for (int idx : pixels) {
      const int i=idx%si, j=idx/si;
      std::vector<double> row(np);
      for(int p=0;p<np;++p) row[p]=BasisValue2D(fittedPeaks2D[grp.members[p]], ivals[i], jvals[j]);
      for(int a=0;a<np;++a) for(int b=0;b<=a;++b) ata[a*np+b]+=row[a]*row[b];
    }
    for(int a=0;a<np;++a) for(int b=0;b<a;++b) ata[b*np+a]=ata[a*np+b];

    // Diagonal of inverse normal matrix, used for conventional LS amplitude ESDs.
    std::vector<double> invdiag(np, 0.0);
    for(int q=0;q<np;++q) {
      std::vector<double> rhs(np,0.0), sol; rhs[q]=1.0;
      if(SolveLinearSystem(ata,rhs,sol,np) && q<(int)sol.size()) invdiag[q]=std::max(0.0,sol[q]);
    }

    for(int z=0; z<nslices; ++z) {
      std::vector<double> aty(np,0.0), coeff;
      for(int idx:pixels) {
        const int i=idx%si, j=idx/si;
        const double yy=stack[z*size+idx];
        for(int p=0;p<np;++p) aty[p]+=BasisValue2D(fittedPeaks2D[grp.members[p]],ivals[i],jvals[j])*yy;
      }
      if(!SolveLinearSystem(ata,aty,coeff,np)) coeff.assign(np,0.0);
      double sse=0.0;
      for(int idx:pixels) {
        const int i=idx%si, j=idx/si; double pred=0.0;
        for(int p=0;p<np;++p) pred+=coeff[p]*BasisValue2D(fittedPeaks2D[grp.members[p]],ivals[i],jvals[j]);
        const double r=stack[z*size+idx]-pred; sse+=r*r;
      }
      const int dof=std::max(1,(int)pixels.size()-np);
      const double var=sse/dof;
      for(int p=0;p<np;++p) {
        const int pm=grp.members[p];
        intensities[pm][z]=coeff[p];
        intensityEsd[pm][z]=std::sqrt(std::max(0.0,var*invdiag[p]));
      }
    }
  }

  if (::mkdir(fitDir.c_str(), 0775) != 0 && errno != EEXIST) return false;
  for(size_t p=0;p<fittedPeaks2D.size();++p) {
    const FitPeak2DLocal& pk=fittedPeaks2D[p];
    const std::string base=fitDir+"/"+(pk.name.empty()?std::to_string(p+1):pk.name);
    FILE* out=std::fopen((base+".out").c_str(),"w");
    if(out) {
      std::fprintf(out,"#\n# %15s%15s\n","Peak Name","Overlap_group");
      std::fprintf(out,"# %15s%15d\n",pk.name.c_str(),pk.group);
      std::fprintf(out,"#\n# --------- Results of the fit -------------\n#\n");
      std::fprintf(out,"# Parameter        Value           Esd      \n");
      const double shapeEsd=std::numeric_limits<double>::quiet_NaN();
      std::fprintf(out,"# %-12s %14.7e %14.7e\n","f01(ppm)",pk.x,shapeEsd);
      std::fprintf(out,"# %-12s %14.7e %14.7e\n","w1(Hz)",pk.sig1*obsI,shapeEsd);
      std::fprintf(out,"# %-12s %14.7e %14.7e\n","g1",pk.voigt1,shapeEsd);
      std::fprintf(out,"# %-12s %14.7e %14.7e\n","f02(ppm)",pk.y,shapeEsd);
      std::fprintf(out,"# %-12s %14.7e %14.7e\n","w2(Hz)",pk.sig2*obsJ,shapeEsd);
      std::fprintf(out,"# %-12s %14.7e %14.7e\n","g2",pk.voigt2,shapeEsd);
      for(int q=0;q<44;++q) std::fputc('#',out); std::fputc('\n',out);
      std::fprintf(out,"# %10s        Intensity      Esd(Int.)\n",zlabel.c_str());
      for(int z=0;z<nslices;++z) std::fprintf(out,"%12.3e   %14.7e %14.7e\n",z<(int)zvals.size()?zvals[z]:z+1.0,intensities[p][z],intensityEsd[p][z]);
      std::fclose(out);
    }
    FILE* dat=std::fopen((base+".dat").c_str(),"w");
    if(dat) {
      std::fprintf(dat,"%11s %11s %11s %11s\n","# F2(ppm)","F1(ppm)","Data","Calc");
      // Match the ordinary 2D FIT .dat writer: use the complete fitted
      // overlap-group ROI rather than testing the single peak ellipse here.
      // BuildEllipsePixels2D also guarantees the same pixel selection/order
      // used by the fit itself.  Emit one blank line between F1 rows and two
      // blank lines between pseudo-Z slices (FUDA/decon convention).
      const int groupIndex = pk.group;
      if (groupIndex < 0 || groupIndex >= static_cast<int>(groups.size()))
      {
        std::cerr << "Protocol3P: invalid stored group " << groupIndex
                  << " for peak " << pk.name << std::endl;
        std::fclose(dat);
        return false;
      }
      const FitGroup2DLocal& outgrp = groups[groupIndex];
      std::vector<int> outPixels;
      BuildEllipsePixels2D(*this, outgrp, fittedPeaks2D, radIppm, radJppm, outPixels);

      for(int z=0;z<nslices;++z) {
        int lastJ = -999999;
        for(size_t qq=0; qq<outPixels.size(); ++qq) {
          const int idx=outPixels[qq];
          const int i=idx%si, j=idx/si;
          if(j != lastJ) {
            if(lastJ != -999999) std::fputc('\n',dat);
            lastJ=j;
          }
          double calc=0.0;
          for(size_t q=0;q<outgrp.members.size();++q) {
            const int pm=outgrp.members[q];
            calc+=intensities[pm][z]*BasisValue2D(fittedPeaks2D[pm],ivals[i],jvals[j]);
          }
          std::fprintf(dat,"%11.4e %11.4e %11.4e %11.4e\n",
                       jvals[j],ivals[i],stack[z*size+idx],calc);
        }
        // Two carriage/blank lines delimit individual pseudo-Z slices.
        std::fputc('\n',dat);
        std::fputc('\n',dat);
      }
      std::fclose(dat);
    }
  }
  return true;
}
