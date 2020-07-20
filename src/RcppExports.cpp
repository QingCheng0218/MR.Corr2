// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "pdsoft.hpp"
#include "GibbsAlpGamEta_ptr.hpp"
#include "ReadGeneFile.hpp"
#include "CalCorr.hpp"
#include "data_loader.hpp"
#include <ctime>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// load_block_file
vector<umat> load_block_file(string block_file);
RcppExport SEXP _MR_Corr2_load_block_file(SEXP block_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< string >::type block_file(block_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(load_block_file(block_file));
    return rcpp_result_gen;
END_RCPP
}
// test_blocks
List test_blocks(arma::ivec bp, arma::ivec chr, std::string block_file);
RcppExport SEXP _MR_Corr2_test_blocks(SEXP bpSEXP, SEXP chrSEXP, SEXP block_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec >::type bp(bpSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< std::string >::type block_file(block_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(test_blocks(bp, chr, block_file));
    return rcpp_result_gen;
END_RCPP
}
// std_setdiff
arma::ivec std_setdiff(arma::ivec& x, arma::ivec& y);
RcppExport SEXP _MR_Corr2_std_setdiff(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(std_setdiff(x, y));
    return rcpp_result_gen;
END_RCPP
}
// LDclump
ivec LDclump(arma::mat& R, double ld_r2_thresh);
RcppExport SEXP _MR_Corr2_LDclump(SEXP RSEXP, SEXP ld_r2_threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type ld_r2_thresh(ld_r2_threshSEXP);
    rcpp_result_gen = Rcpp::wrap(LDclump(R, ld_r2_thresh));
    return rcpp_result_gen;
END_RCPP
}
// Cal_blockR
List Cal_blockR(arma::ivec& bp, arma::ivec& chr, arma::uvec& avbIndex, arma::uvec& idx4panel, std::string block_file, std::string stringname3, double ld_r2_thresh, int coreNum, double lam);
RcppExport SEXP _MR_Corr2_Cal_blockR(SEXP bpSEXP, SEXP chrSEXP, SEXP avbIndexSEXP, SEXP idx4panelSEXP, SEXP block_fileSEXP, SEXP stringname3SEXP, SEXP ld_r2_threshSEXP, SEXP coreNumSEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type bp(bpSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type avbIndex(avbIndexSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type idx4panel(idx4panelSEXP);
    Rcpp::traits::input_parameter< std::string >::type block_file(block_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type stringname3(stringname3SEXP);
    Rcpp::traits::input_parameter< double >::type ld_r2_thresh(ld_r2_threshSEXP);
    Rcpp::traits::input_parameter< int >::type coreNum(coreNumSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, ld_r2_thresh, coreNum, lam));
    return rcpp_result_gen;
END_RCPP
}
// Cal_block_Rmatrix
List Cal_block_Rmatrix(arma::ivec& bp, arma::ivec& chr, arma::uvec& avbIndex, arma::uvec& idx4panel, std::string block_file, std::string stringname3, double ld_r2_thresh, int coreNum, double lam);
RcppExport SEXP _MR_Corr2_Cal_block_Rmatrix(SEXP bpSEXP, SEXP chrSEXP, SEXP avbIndexSEXP, SEXP idx4panelSEXP, SEXP block_fileSEXP, SEXP stringname3SEXP, SEXP ld_r2_threshSEXP, SEXP coreNumSEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type bp(bpSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type avbIndex(avbIndexSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type idx4panel(idx4panelSEXP);
    Rcpp::traits::input_parameter< std::string >::type block_file(block_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type stringname3(stringname3SEXP);
    Rcpp::traits::input_parameter< double >::type ld_r2_thresh(ld_r2_threshSEXP);
    Rcpp::traits::input_parameter< int >::type coreNum(coreNumSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(Cal_block_Rmatrix(bp, chr, avbIndex, idx4panel, block_file, stringname3, ld_r2_thresh, coreNum, lam));
    return rcpp_result_gen;
END_RCPP
}
// Cal_blockinf
List Cal_blockinf(arma::ivec& bp, arma::ivec& chr, std::string block_file);
RcppExport SEXP _MR_Corr2_Cal_blockinf(SEXP bpSEXP, SEXP chrSEXP, SEXP block_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type bp(bpSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< std::string >::type block_file(block_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(Cal_blockinf(bp, chr, block_file));
    return rcpp_result_gen;
END_RCPP
}
// MRcorr
List MRcorr(arma::vec& gammah, arma::vec& Gammah, arma::vec& se1, arma::vec& se2, SEXP opts);
RcppExport SEXP _MR_Corr2_MRcorr(SEXP gammahSEXP, SEXP GammahSEXP, SEXP se1SEXP, SEXP se2SEXP, SEXP optsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type gammah(gammahSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Gammah(GammahSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se1(se1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2(se2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
    rcpp_result_gen = Rcpp::wrap(MRcorr(gammah, Gammah, se1, se2, opts));
    return rcpp_result_gen;
END_RCPP
}
// MRCorr2Sim
List MRCorr2Sim(arma::vec& gammah, arma::vec& Gammah, arma::vec& se1, arma::vec& se2, arma::mat R, arma::umat block_inf, int coreNum, SEXP opts);
RcppExport SEXP _MR_Corr2_MRCorr2Sim(SEXP gammahSEXP, SEXP GammahSEXP, SEXP se1SEXP, SEXP se2SEXP, SEXP RSEXP, SEXP block_infSEXP, SEXP coreNumSEXP, SEXP optsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type gammah(gammahSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Gammah(GammahSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se1(se1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2(se2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type block_inf(block_infSEXP);
    Rcpp::traits::input_parameter< int >::type coreNum(coreNumSEXP);
    Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
    rcpp_result_gen = Rcpp::wrap(MRCorr2Sim(gammah, Gammah, se1, se2, R, block_inf, coreNum, opts));
    return rcpp_result_gen;
END_RCPP
}
// MRCorr2Real
Rcpp::List MRCorr2Real(arma::ivec& bp, arma::ivec& chr, arma::uvec& avbIndex, arma::uvec& idx4panel, std::string& block_file, std::string stringname3, double ld_r2_thresh, arma::vec& bh1, arma::vec& bh2, arma::vec& se1, arma::vec& se2, double lam, int coreNum, SEXP opts);
RcppExport SEXP _MR_Corr2_MRCorr2Real(SEXP bpSEXP, SEXP chrSEXP, SEXP avbIndexSEXP, SEXP idx4panelSEXP, SEXP block_fileSEXP, SEXP stringname3SEXP, SEXP ld_r2_threshSEXP, SEXP bh1SEXP, SEXP bh2SEXP, SEXP se1SEXP, SEXP se2SEXP, SEXP lamSEXP, SEXP coreNumSEXP, SEXP optsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type bp(bpSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type avbIndex(avbIndexSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type idx4panel(idx4panelSEXP);
    Rcpp::traits::input_parameter< std::string& >::type block_file(block_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type stringname3(stringname3SEXP);
    Rcpp::traits::input_parameter< double >::type ld_r2_thresh(ld_r2_threshSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type bh1(bh1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type bh2(bh2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se1(se1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2(se2SEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< int >::type coreNum(coreNumSEXP);
    Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
    rcpp_result_gen = Rcpp::wrap(MRCorr2Real(bp, chr, avbIndex, idx4panel, block_file, stringname3, ld_r2_thresh, bh1, bh2, se1, se2, lam, coreNum, opts));
    return rcpp_result_gen;
END_RCPP
}
// Cal_block_SimR
mat Cal_block_SimR(umat block_inf, arma::umat& X, double lam);
RcppExport SEXP _MR_Corr2_Cal_block_SimR(SEXP block_infSEXP, SEXP XSEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< umat >::type block_inf(block_infSEXP);
    Rcpp::traits::input_parameter< arma::umat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(Cal_block_SimR(block_inf, X, lam));
    return rcpp_result_gen;
END_RCPP
}
// fastSigLm
List fastSigLm(const arma::vec& y, const arma::mat& X);
RcppExport SEXP _MR_Corr2_fastSigLm(SEXP ySEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastSigLm(y, X));
    return rcpp_result_gen;
END_RCPP
}
// getLineNum
int getLineNum(std::string filename);
RcppExport SEXP _MR_Corr2_getLineNum(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(getLineNum(filename));
    return rcpp_result_gen;
END_RCPP
}
// ReadSNPinfo
Rcpp::List ReadSNPinfo(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname, IntegerVector chr, IntegerVector bp, NumericVector morgan, int N);
RcppExport SEXP _MR_Corr2_ReadSNPinfo(SEXP stringnameSEXP, SEXP A1SEXP, SEXP A2SEXP, SEXP rsnameSEXP, SEXP chrSEXP, SEXP bpSEXP, SEXP morganSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type stringname(stringnameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type A2(A2SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type rsname(rsnameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type bp(bpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadSNPinfo(stringname, A1, A2, rsname, chr, bp, morgan, N));
    return rcpp_result_gen;
END_RCPP
}
// Read_summarystat
void Read_summarystat(std::string stringname, IntegerVector SA1, IntegerVector SA2, CharacterVector rsname, NumericVector betah, NumericVector s2, NumericVector pvalue, IntegerVector chr, IntegerVector bp, int N);
RcppExport SEXP _MR_Corr2_Read_summarystat(SEXP stringnameSEXP, SEXP SA1SEXP, SEXP SA2SEXP, SEXP rsnameSEXP, SEXP betahSEXP, SEXP s2SEXP, SEXP pvalueSEXP, SEXP chrSEXP, SEXP bpSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type stringname(stringnameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type SA1(SA1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type SA2(SA2SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type rsname(rsnameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betah(betahSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pvalue(pvalueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type bp(bpSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Read_summarystat(stringname, SA1, SA2, rsname, betah, s2, pvalue, chr, bp, N);
    return R_NilValue;
END_RCPP
}
// select
CharacterVector select(CharacterVector vec_, NumericVector idx_);
RcppExport SEXP _MR_Corr2_select(SEXP vec_SEXP, SEXP idx_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type vec_(vec_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type idx_(idx_SEXP);
    rcpp_result_gen = Rcpp::wrap(select(vec_, idx_));
    return rcpp_result_gen;
END_RCPP
}
// matchsnp
Rcpp::List matchsnp(std::string stringname1, std::string stringname2, std::string stringname3, bool matchExp);
RcppExport SEXP _MR_Corr2_matchsnp(SEXP stringname1SEXP, SEXP stringname2SEXP, SEXP stringname3SEXP, SEXP matchExpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type stringname1(stringname1SEXP);
    Rcpp::traits::input_parameter< std::string >::type stringname2(stringname2SEXP);
    Rcpp::traits::input_parameter< std::string >::type stringname3(stringname3SEXP);
    Rcpp::traits::input_parameter< bool >::type matchExp(matchExpSEXP);
    rcpp_result_gen = Rcpp::wrap(matchsnp(stringname1, stringname2, stringname3, matchExp));
    return rcpp_result_gen;
END_RCPP
}
// matchscreen
Rcpp::List matchscreen(std::string screenname, std::string stringname1, std::string stringname2, std::string stringname3, double pva_cutoff, bool matchExp);
RcppExport SEXP _MR_Corr2_matchscreen(SEXP screennameSEXP, SEXP stringname1SEXP, SEXP stringname2SEXP, SEXP stringname3SEXP, SEXP pva_cutoffSEXP, SEXP matchExpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type screenname(screennameSEXP);
    Rcpp::traits::input_parameter< std::string >::type stringname1(stringname1SEXP);
    Rcpp::traits::input_parameter< std::string >::type stringname2(stringname2SEXP);
    Rcpp::traits::input_parameter< std::string >::type stringname3(stringname3SEXP);
    Rcpp::traits::input_parameter< double >::type pva_cutoff(pva_cutoffSEXP);
    Rcpp::traits::input_parameter< bool >::type matchExp(matchExpSEXP);
    rcpp_result_gen = Rcpp::wrap(matchscreen(screenname, stringname1, stringname2, stringname3, pva_cutoff, matchExp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MR_Corr2_load_block_file", (DL_FUNC) &_MR_Corr2_load_block_file, 1},
    {"_MR_Corr2_test_blocks", (DL_FUNC) &_MR_Corr2_test_blocks, 3},
    {"_MR_Corr2_std_setdiff", (DL_FUNC) &_MR_Corr2_std_setdiff, 2},
    {"_MR_Corr2_LDclump", (DL_FUNC) &_MR_Corr2_LDclump, 2},
    {"_MR_Corr2_Cal_blockR", (DL_FUNC) &_MR_Corr2_Cal_blockR, 9},
    {"_MR_Corr2_Cal_block_Rmatrix", (DL_FUNC) &_MR_Corr2_Cal_block_Rmatrix, 9},
    {"_MR_Corr2_Cal_blockinf", (DL_FUNC) &_MR_Corr2_Cal_blockinf, 3},
    {"_MR_Corr2_MRcorr", (DL_FUNC) &_MR_Corr2_MRcorr, 5},
    {"_MR_Corr2_MRCorr2Sim", (DL_FUNC) &_MR_Corr2_MRCorr2Sim, 8},
    {"_MR_Corr2_MRCorr2Real", (DL_FUNC) &_MR_Corr2_MRCorr2Real, 14},
    {"_MR_Corr2_Cal_block_SimR", (DL_FUNC) &_MR_Corr2_Cal_block_SimR, 3},
    {"_MR_Corr2_fastSigLm", (DL_FUNC) &_MR_Corr2_fastSigLm, 2},
    {"_MR_Corr2_getLineNum", (DL_FUNC) &_MR_Corr2_getLineNum, 1},
    {"_MR_Corr2_ReadSNPinfo", (DL_FUNC) &_MR_Corr2_ReadSNPinfo, 8},
    {"_MR_Corr2_Read_summarystat", (DL_FUNC) &_MR_Corr2_Read_summarystat, 10},
    {"_MR_Corr2_select", (DL_FUNC) &_MR_Corr2_select, 2},
    {"_MR_Corr2_matchsnp", (DL_FUNC) &_MR_Corr2_matchsnp, 4},
    {"_MR_Corr2_matchscreen", (DL_FUNC) &_MR_Corr2_matchscreen, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_MR_Corr2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
