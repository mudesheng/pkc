package discretecosinetransform

// #cgo CFLAGS: -I/usr/local/include -I/usr/include/lapacke/ -std=c99
// #cgo LDFLAGS: -L ./ -lm -lcv -lcxcore -ljpeg -lpng -lhighgui -llapack -llapacke
// #include <cv.h>
// #include <stdio.h>
// #include <math.h>
// #include <highgui.h>
// #include <cxcore.h>
// #include <lapacke.h>
/*


typedef struct {
  double *s;
  double *u;
  double *vt;
} SV ;

typedef struct {
  double *cdctarr;
  CvMat *dst;
} DCTRet;

void jpgcvDCT() {
  IplImage *input, *output, *b;
  input = cvLoadImage("dct_orig.png",CV_LOAD_IMAGE_COLOR);
  int width = input->width;
  int height = input->height;
  // b = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U,1);
  b = cvCreateImage(cvSize(width, height), CV_32FC1,1);
  IplImage *b_dct;
  b_dct=cvCreateImage(cvSize(width,height),CV_32FC1,1);
  cvDCT(b,b_dct,0);  // doubt??
  cvShowImage("blue",b_dct);
  cvSaveImage("/db/goworkspace/src/pkc/dct.jpg", b_dct);
}

double* NewdoubleArray(int length) {
  double *arr = calloc(length, sizeof(double));
  return arr;
}

inline void SetdoubleArray(double *arr, int pos, double val) {
  arr[pos] = val;
}

// Auxiliary routine: printing a matrix
void print_matrix( char* desc, int m, int n, double* a, int lda  ) {
	int i, j;
	printf( "\n %s\n", desc  );
	for( i = 0; i < m; i++  ) {
		for( j = 0; j < n; j++  ) printf( " %6.2f", a[i*lda+j]);
		printf( "\n"  );
	}
}

void print_cvMat(CvMat *mat) {
  printf("[print_cvMat]cvMat:\n");
  for(int i = 0; i <mat->rows;i++) {
    for(int j = 0; j < mat->cols; j++) {
      printf(" %6.2f", cvmGet(mat, i, j));
    }
    printf("\n");
  }
}

DCTRet ReadsDCT(double *data, const int length) {
  // printf("data pointer address: %p\t", data);
  // printf("length: %d, clang arr:\t", length);
  // printf("%d:%f, %d:%f\n", 0, data[0], length-1, data[length-1]);

  // for( int i = 0; i < length; i++) {
  //   printf("%f ");
  // }
  // printf("\n");
  // double row = sqrt(length);
  // if(row * row != length  && (int)row % 2 == 0) {
  //   fprintf(stderr, "[ReadsDCT]the input data matrix rows and columns must be equaled and even, exit...\n");
  //     exit(1);
  // }

  CvMat src64 ;
  CvMat *dst64 = cvCreateMat(1, length, CV_64FC1);
  cvInitMatHeader(&src64, 1, length ,CV_64FC1, data, CV_AUTOSTEP);
  cvDCT(&src64,dst64,CV_DXT_FORWARD | CV_DXT_ROWS);
  // CvMat *dst64 = cvCreateMat(row, row, CV_64FC1);
  // cvInitMatHeader(&src64, row, row ,CV_64FC1, data, CV_AUTOSTEP);
  // cvDCT(&src64,dst64,CV_DXT_FORWARD);
  // print_cvMat(dst64);
  // print_matrix("dst64", dst64->rows, dst64->cols, (double*)dst64->data.db, dst64->cols);
  // free work
  cvReleaseData(&src64);
  // free(data);

  DCTRet dctRet;
  dctRet.cdctarr = (double*)dst64->data.db;
  dctRet.dst = dst64;
  return dctRet;

}

// SV PSVD() {
//   pdgesvd('V','N',m, n, dataMat, 1,1,)
// }

SV SVD(double *dataMat, const int row, const int column) {
	// Locals
	int M = row, N = column;
	int LDA = N, LDU = M, LDVT = N;
	int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT,info;

	int min = row > column ? column :row;
	// double superb[min-1];
	double *superb = calloc(min-1, sizeof(double));
	// double s[N], u[LDU *M], vt[LDVT *N];
	double *s = calloc(N, sizeof(double));
	// double *u = calloc(LDU * M, sizeof(double));
	double *u = calloc(LDU * M, sizeof(double));
	double *vt = calloc(LDVT *N, sizeof(double));

	// Executable statements
	printf( "[C.SVD] DGESVD Example Program Results\n");
	fflush(stdout);

	// Compute SVD
	// info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR,'A', 'A', m, n, dataMat, lda, s, u, ldu, vt, ldvt,superb);
	info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR,'O', 'N', m, n, dataMat, lda, s, u, ldu, NULL, ldvt,superb);
	// dgesvd("All", "All", &m, &n, dataMat, &lda, s, u, &ldu, vt, &ldvt, wkopt, &lwork);
	// lwork = (int)wkopt;
	// work = (double*)malloc(lwork * sizeof(double));
	// LAPACKE_sgesvd("All", "All", &m, &n, dataMat, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
	// // check for convergence
	if(info > 0) {
	  fprintf(stderr, "[__FUNC__] The algorithm computing SVD failed to converge.\n");
	  fflush(stderr);
	  exit(1);
	}
	// Print singular values
	print_matrix( "Singular values", 1, n, s, 1);
	// Print left singular vectors
	print_matrix( "Left singular vectors (stored columnwise)", m, n, u, ldu  );
	// Print right singular vectors
	print_matrix( "Right singular vectors (stored rowwise)", n, n, vt, ldvt  );
	// Free workspace
	free(superb);
	// free(s); free(u); free(vt);
	SV sv;
	sv.s = s; sv.u = u; sv.vt = vt;

	return sv;
}

  // IplImage *image=0;
  // IplImage *image1=0;
  //  IplImage *image2=0;
  // IplImage *image3=0;
  // CvMat* src32 = NULL;
  // CvMat arr;
  //       CvMat* dst=NULL;
  //       CvMat* idst=NULL;
	// double data[] = {1,2,3,4,\
			    // 6,7,8,10,\
			    // 11,12,14,15,\
			    // 16,17,19,20,\
			    // 21,22,24,25};

	// image=cvLoadImage("timeline.jpg",CV_LOAD_IMAGE_COLOR);

	  // image1=cvCreateImage(cvSize(image->width,image->height),IPL_DEPTH_32F,1);
	    // image2=cvCreateImage(cvSize(image->width, image->height),IPL_DEPTH_32F,1);
	      // image3=cvCreateImage(cvSize(image->width, image->height),IPL_DEPTH_32F,1);
	      // src32=cvCreateMat(image->height, image->width,CV_32FC1);
	      // dst=cvCreateMat(image->height, image->width, CV_64FC1);
	      // idst=cvCreateMat(image->height, image->width, CV_64FC1);

	      // // cvCvtColor(image,image1,CV_YCrCb2BGR);
	      // cvInitMatHeader(&arr, 4,4,CV_64FC1, data, 0);
	      // // cvCvtColor(image,image1,CV_RGB2Luv);
	         // // cvConvertScale(image1, src32,1,0);
		    // cvDCT(&arr,dst,CV_DXT_FORWARD);
		    // // cvDCT(im8,image2,CV_DXT_FORWARD);
		        // // cvConvertImage(image2,dst,0);
			// printf("dst:\n");
			// int i;
			// int h =4, w = 4;
			// for(i=0;i<h;i++)
			// {
			  // int j;
			   // for(j=0;j<w;j++)
			   // {
			       // printf("%5d",abs(cvmGet(dst,i,j)));

			   // }
			    // printf("\n");

			// }
			// cvDCT(dst,idst,CV_DXT_INVERSE);
			// printf("idst:\n");
			// for(i=0;i<h;i++)
			// {
			  // int j;
			   // for(j=0;j<w;j++)
			   // {
			       // printf("%5d",abs(cvmGet(idst,i,j)));

			   // }
			    // printf("\n");

			// }
			// // cvNamedWindow("jpeg",1);
			// //    cvNamedWindow("gray",1);
			// //       cvNamedWindow("dct",1);
			// //       cvNamedWindow("idct",1);
			          // // cvShowImage("jpeg",image);
				     // // cvShowImage("gray",image1);
				        // // cvShowImage("dct",image2);
					   // // cvShowImage("idct",image3);
					   // // cvSaveImage("dct.jpg", idst);
					       // // cvWaitKey(0);
					           // cvReleaseImage(&image);
						   // cvReleaseImage(&image1);
						   // cvReleaseImage(&image2);
						   // cvReleaseImage(&image3);
						   // cvReleaseMat(&src32);
						   // cvReleaseMat(&dst);

// }
*/
import "C"

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"os/exec"
	"sort"
	"strconv"
	"strings"
	"time"
	"unsafe"

	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
	// smatrix "github.com/skelterjohn/go.matrix"
	// "github.com/gonum/matrix/mat64"
)

type DCTArr struct {
	dct []float64
}

type SvdVect struct {
	ID       string
	Position int
	LeftSV   []float64
}

type DigitSrc struct {
	ID  string
	src []C.double
}

type ReadIDCount struct {
	ID    string
	count int
}

type ReadIDCountArr []ReadIDCount

type NearestReadInfo struct {
	Pos       int    // reference read start position
	ID1       string // prime nearest read ID
	StartPos1 int    // prime read start position
	ID2       string // second nearest read ID
	StartPos2 int    // second nearest read start position
}

var totalDeviation float64 = 0.02
var singleDeviation float64 = 0.1

var TopSigmaArr []float64

var Depth int = 20
var SegmentLen int = 8192
var StepLen int = 100

type SvdVectArr []SvdVect

func (sva SvdVectArr) Len() int {
	return len(sva)
}

func (sva SvdVectArr) Less(i, j int) bool {
	isum := math.Abs(sva[i].LeftSV[0] * TopSigmaArr[0])
	jsum := math.Abs(sva[j].LeftSV[0] * TopSigmaArr[0])
	for x := 0; x < len(sva[i].LeftSV)-1 && x < len(sva[j].LeftSV)-1; x += 1 {
		a, b := sva[i].LeftSV[x], sva[j].LeftSV[x]
		if math.Abs((a-b)/((a+b)/2)) < singleDeviation &&
			(a*b > 0 && sva[i].LeftSV[x+1]*sva[j].LeftSV[x+1] > 0) &&
			math.Abs(a-b)*TopSigmaArr[x] <
				math.Abs(sva[i].LeftSV[x+1]-sva[j].LeftSV[x+1])*TopSigmaArr[x+1] {

			isum += math.Abs(sva[i].LeftSV[x+1] * TopSigmaArr[x+1])
			jsum += math.Abs(sva[i].LeftSV[x+1] * TopSigmaArr[x+1])
			if math.Abs(isum-jsum) < TopSigmaArr[0]*totalDeviation {
				continue
			} else {
				return isum < jsum
			}
		} else {
			return sva[i].LeftSV[x] < sva[j].LeftSV[x]
		}
	}

	return len(sva[i].LeftSV) < len(sva[j].LeftSV)
}

func (sva SvdVectArr) Swap(i, j int) {
	sva[i], sva[j] = sva[j], sva[i]
}

// var DNA2Num = [5]C.double{0, 63, 127, 191, 255}

var DNA2Num = [5]C.double{0, -1.5, -0.5, 0.5, 1.5}

func DNA2Pixel(tl alphabet.Letters) (srcarr []C.double) {
	// the size DCT need must be [even*even] number Matrix
	// var size int
	// sq := math.Sqrt(float64(len(tl)))
	// m := int(math.Floor(sq))
	// if m%2 != 0 {
	// 	m -= 1
	// }

	// srcarr = make([]C.double, len(tl))
	// size = m * m

	// size := min
	// if len(tl) < min {
	// 	if len(tl)%2 != 0 {
	// 		size = len(tl) - 1
	// 	} else {
	// 		size = len(tl)
	// 	}
	// }

	// stl := tl.String()
	for _, e := range tl {
		// fmt.Printf("e:%v\n", e)
		switch e {
		case 'A':
			srcarr = append(srcarr, DNA2Num[1])
		case 'C':
			srcarr = append(srcarr, DNA2Num[2])
		case 'G':
			srcarr = append(srcarr, DNA2Num[3])
		case 'T':
			srcarr = append(srcarr, DNA2Num[4])
		default:
			srcarr = append(srcarr, DNA2Num[0])
		}
	}

	// fmt.Printf("tl: %v\nsrcarr:%v\n", tl.String(), srcarr)

	return
}

func DNA2PixelRev(tl alphabet.Letters) (srcarr []C.double) {
	for i := len(tl) - 1; i >= 0; i -= 1 {
		// fmt.Printf("e:%v\n", e)
		switch tl[i] {
		case 'A':
			srcarr = append(srcarr, DNA2Num[4])
		case 'C':
			srcarr = append(srcarr, DNA2Num[3])
		case 'G':
			srcarr = append(srcarr, DNA2Num[2])
		case 'T':
			srcarr = append(srcarr, DNA2Num[1])
		default:
			srcarr = append(srcarr, DNA2Num[0])
		}
	}

	// fmt.Printf("tl: %v\nsrcarr:%v\n", tl.String(), srcarr)

	return
}

func AdjustSrcArr(dctarr []C.double) []C.double {
	sq := math.Sqrt(float64(len(dctarr)))
	if sq-math.Floor(sq) != 0 || int(sq)%2 != 0 {
		size := int(math.Floor(sq)) + 1
		if size%2 != 0 {
			size += 1
		}
		size *= size
		for i := len(dctarr); i < size; i++ {
			dctarr = append(dctarr, DNA2Num[0])
		}
	}

	return dctarr
}

func Gofloat64Arr2CdoubleArr(srcarr []float64, csrcarr unsafe.Pointer) {
	p := uintptr(csrcarr)
	for _, e := range srcarr {
		// *(*C.double)(unsafe.Pointer(&p)) = C.double(e)
		*(*C.double)(unsafe.Pointer(&p)) = C.double(e)
		fmt.Printf("%v\t", p)
		fmt.Printf("%f\t", *(*float64)(csrcarr))
		p += unsafe.Sizeof(e)
		fmt.Printf("%v\t", p)
	}
	fmt.Printf("%f\t", *(*float64)(csrcarr))
}

func CdoubleArr2Gofloat64Arr(cdctarr unsafe.Pointer, size int) (dctarr []float64) {
	p := uintptr(cdctarr)

	for i := 0; i < size; i++ {
		j := *(*float64)(unsafe.Pointer(p))
		dctarr = append(dctarr, j)
		p += unsafe.Sizeof(j)
	}

	return dctarr
}

func PrintCdctarr(cdctarr *C.double, length int) {
	p := uintptr(unsafe.Pointer(cdctarr))
	row := int(math.Sqrt(float64(length)))
	for i := 0; i < row; i++ {
		for j := 0; j < row; j++ {
			x := *(*float64)(unsafe.Pointer(p))
			fmt.Printf("%f ", x)
			p += unsafe.Sizeof(x)
		}
		fmt.Println()
	}
}

// get the most high power(left-up) square matrix
func GetPrimeDCTArr(dctarr []float64, length int, rows, cols int) (primedct []float64) {
	if rows*cols >= length {
		log.Fatalf("[GetPrimeDCTArr] Get size: %dX%d matrix bigger than length:%d\n", rows, cols, length)
	}

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			primedct = append(primedct, dctarr[i*cols+j])
		}
	}

	return primedct
}

func CheckSuccAndListDir(dirPath string) (files []string, err error) {
	if _, err := os.Stat(dirPath + string(os.PathSeparator) + "_SUCCESS"); err != nil {
		return nil, err
	}

	dir, err := ioutil.ReadDir(dirPath)
	if err != nil {
		return nil, err
	}

	for _, fi := range dir {
		if fi.IsDir() {
			continue
		}
		if strings.HasPrefix(fi.Name(), "part") {
			files = append(files, dirPath+string(os.PathSeparator)+fi.Name())
		}
	}

	return files, nil
}

func StrArr2float64Arr(inArr []string) (outArr []float64) {
	for _, e := range inArr {
		fe, _ := strconv.ParseFloat(e, 64)
		outArr = append(outArr, fe)
	}

	return outArr
}

func PrintPrimeDCT(svdVect SvdVect, primedct []float64, rows, cols int) {
	fmt.Printf(">%s/%d\n", svdVect.ID, svdVect.Position)
	fmt.Fprintf(os.Stderr, "['%d', %6.2f],\n", svdVect.Position, primedct[0])
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			fmt.Printf("%.0f ", primedct[i*cols+j])
		}
		fmt.Printf("\n")
	}
}

func CompareDctArr(lastDctArr, primedct []float64) {
	var difNum int = 0
	var difsign int = 0
	for i := 0; i < len(lastDctArr); i++ {
		v1, v2 := lastDctArr[i], primedct[i]
		if lastDctArr[i]*primedct[i] < 0 {
			v1 = math.Abs(v1)
			v2 = math.Abs(v2)
			difsign++
		}
		if math.Abs(v1-v2) > 1 {
			difNum++
		}
	}

	fmt.Printf("difNum: %d/%d, difSign percent: %d\n", difNum, len(lastDctArr), difsign*100/len(lastDctArr))
}

// func ReadIDCountArrIndex(rIDCarr []ReadIDCount, valk string) int {
func (rIDCarr ReadIDCountArr) Index(valk string) int {
	for i, item := range rIDCarr {
		if item.ID == valk {
			return i
		}
	}

	return -1
}

func InsertReadIDCount(IDcount map[string]ReadIDCountArr, srck, valk string) {
	var rIDC ReadIDCount
	rIDC.ID = valk
	rIDC.count = 1
	if arr, ok := IDcount[srck]; ok {
		pos := arr.Index(valk)
		if pos >= 0 {
			IDcount[srck][pos].count += 1
		} else {
			IDcount[srck] = append(arr, rIDC)
		}
	} else {
		var rIDCarr ReadIDCountArr
		IDcount[srck] = append(rIDCarr, rIDC)
	}
}

func CleanReadIDCountArr(v ReadIDCountArr, MinCount int) (nv ReadIDCountArr) {
	for _, item := range v {
		if item.count >= MinCount {
			nv = append(nv, item)
		}
	}

	return nv
}

func (rca ReadIDCountArr) Len() int {
	return len(rca)
}

func (rca ReadIDCountArr) Swap(i, j int) {
	rca[i], rca[j] = rca[j], rca[i]
}

// decrease sort
func (rca ReadIDCountArr) Less(i, j int) bool {
	return rca[i].count > rca[j].count
}

func ConstrctIDcountFromfn(sortedIDRelationfn string, MinCount int) map[string]ReadIDCountArr {
	IDcount := make(map[string]ReadIDCountArr)
	sortedIDRelationfp, err := os.Open(sortedIDRelationfn)
	defer sortedIDRelationfp.Close()
	if err != nil {
		log.Fatalf("[ConstrctIDcountFromfn] open %s failed, err: %v\n", sortedIDRelationfn, err)
	}
	buffp := bufio.NewReader(sortedIDRelationfp)

	eof := false
	for !eof {
		line, err := buffp.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				err = nil
				eof = true
				continue
			} else {
				log.Fatalf("[ConstrctIDcountFromfn] read file: %s failed, err : %v\n", sortedIDRelationfn, err)
			}
		}

		fds := strings.Fields(line)
		count, err := strconv.Atoi(fds[0])
		if err != nil {
			log.Fatalf("[ConstrctIDcountFromfn] err: %v\n", err)
		}
		if count >= MinCount {
			var rIDC ReadIDCount
			rIDC.ID = fds[2]
			rIDC.count = count
			IDcount[fds[1]] = append(IDcount[fds[1]], rIDC)
			// rIDC.ID = fds[1]
			// IDcount[fds[2]] = append(IDcount[fds[2]], rIDC)
		}
	}

	return IDcount
}

func IsInStringArr(storedID []string, ID string) bool {
	for _, item := range storedID {
		if item == ID {
			return true
		}
	}

	return false
}

func GetUnionArr(v, rv ReadIDCountArr) (unionArr ReadIDCountArr) {
	vMaxSize := Depth
	if vMaxSize > len(v) {
		vMaxSize = len(v)
	}
	for i, e := range v[:vMaxSize] {
		if e.ID[len(e.ID)-2:] == "RC" {
			v[i].ID = e.ID[:len(e.ID)-3]
		}
	}

	rvMaxSize := Depth
	if rvMaxSize > len(rv) {
		rvMaxSize = len(rv)
	}

	for i, e := range rv[:rvMaxSize] {
		if e.ID[len(e.ID)-2:] == "RC" {
			rv[i].ID = e.ID[:len(e.ID)-3]
		}
	}

	for _, ei := range v[:vMaxSize] {
		for _, ej := range rv[:rvMaxSize] {
			if ei.ID == ej.ID {
				var rIDC ReadIDCount
				rIDC.ID = ei.ID
				rIDC.count = ei.count + ej.count
				unionArr = append(unionArr, rIDC)
				break
			}
		}
	}

	sort.Sort(unionArr)

	return unionArr
}

func ConstructNearestReadRelation(svdVectArr SvdVectArr, readLenMap map[string]int) (NearestReadInfoMap map[string]NearestReadInfo) {
	// extract all near read ID and count occur number
	maxWinSize := 40
	// IDcount := make(map[string]ReadIDCountArr)
	IDRelationfn := "IDrelation.txt"
	IDRelationfp, err := os.Create(IDRelationfn)
	if err != nil {
		log.Fatalf("[ConstructNearestReadRelation] create %s err: %v\n", IDRelationfn, err)
	}

	for i, vect := range svdVectArr {
		// sep_i := strings.LastIndex(vect.ID, "/")
		// fmt.Printf("%s\n", vect.ID)
		IDi := vect.ID
		var rIDi string
		if IDi[len(IDi)-2:] == "RC" {
			rIDi = IDi[:len(IDi)-3]
		} else {
			rIDi = IDi + "_RC"
		}
		leftExtend, rightExtend := true, true
		var storedID []string
		for j := 1; j < maxWinSize; j += 1 {
			if leftExtend == false && rightExtend == false {
				break
			}
			if leftExtend == true && i-j >= 0 {
				vj := svdVectArr[i-j]
				if vj.ID != IDi && vj.ID != rIDi {
					if IsInStringArr(storedID, vj.ID) == false {
						fmt.Fprintf(IDRelationfp, "%s\t%s\t%d\t%d\n", IDi, vj.ID, vect.Position, vj.Position)
						storedID = append(storedID, vj.ID)
					}
				} else {
					leftExtend = false
				}
			}

			if rightExtend == true && i+j < len(svdVectArr) {
				vj := svdVectArr[i+j]
				if vj.ID != IDi && vj.ID != rIDi {
					if IsInStringArr(storedID, vj.ID) == false {
						fmt.Fprintf(IDRelationfp, "%s\t%s\t%d\t%d\n", IDi, vj.ID, vect.Position, vj.Position)
						storedID = append(storedID, vj.ID)
					}
				} else {
					rightExtend = false
				}
			}

		}
		/*
			for j := i - 1; j >= 0 && j > i-maxWinSize; j -= 1 {
				vj := svdVectArr[j]
				// sep_j := strings.LastIndex(vj.ID, "/")
				IDj := vj.ID
				// add ReadIDCount relation to the IDCount
				if IDi != IDj {
					if IDi > IDj {
						IDi, IDj = IDj, IDi
					}
					fmt.Fprintf(IDRelationfp, "%s\t%s\n", IDi, IDj)
					// InsertReadIDCount(IDcount, IDi, IDj)
					// InsertReadIDCount(IDcount, IDj, IDi)
				}
			} */
	}

	IDRelationfp.Close()
	// call system sort and uniq command
	start := time.Now()
	fmt.Printf("Begin sort ReadIDCountArr....")
	// cut the first two fields
	cutfn := IDRelationfn + ".cutf2"
	// cutfcmd := exec.Command("cut", "-f", "1", IDRelationfn, ">", cutfn)
	cutcmd := exec.Command("cut", "-f", "1,2", IDRelationfn)
	cutfp, err := os.Create(cutfn)
	if err != nil {
		log.Fatalf("[ConstructNearestReadRelation] open %s err: %v\n", cutfn, err)
	}
	cutcmd.Stdout = cutfp
	if err := cutcmd.Run(); err != nil {
		log.Fatalf("[ConstructNearestReadRelation] cut Command called error: %v\n", err)
	}
	cutfp.Close()

	sortedIDRelationfn := cutfn + ".sorted"
	sortcmd := exec.Command("sort", "--parallel", "2", "-S", "4G", "-k", "1,1", "-k", "2,2", "-k", "3,3n", "-o", sortedIDRelationfn, IDRelationfn)
	uniqcmd := exec.Command("uniq", "-dc", "-", sortedIDRelationfn)
	reader, writer := io.Pipe()
	sortcmd.Stdout = writer
	uniqcmd.Stdin = reader
	err1 := sortcmd.Start()
	err2 := uniqcmd.Start()
	err3 := sortcmd.Wait()
	writer.Close()
	err4 := uniqcmd.Wait()

	if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
		log.Fatalf("[ConstructNearestReadRelation] system sort command called error: %v, %v, %v, %v\n", err1, err2, err3, err4)
	}
	fmt.Printf("finished, used %vm\n", time.Now().Sub(start).Minutes())

	// Read from sortedIDRelationfn
	MinCount := 5
	IDcount := ConstrctIDcountFromfn(sortedIDRelationfn, MinCount)

	// delete lower than cutoff count and sort ReadIDCountArr
	for _, v := range IDcount {
		// IDcount[k] = CleanReadIDCountArr(v, MinCount)
		sort.Sort(v)
	}

	// sort ReadIDCountArr
	// for _, v := range IDcount {
	// 	sort.Sort(v)
	// }

	// find most probable nearest read
	var symmetryMatch int
	var totalNum int
	var nounionNum int
	for k, v := range IDcount {
		// fmt.Printf("%v\n", v)
		if k[len(k)-2:] != "RC" {
			totalNum += 1
			rk := k + "_RC"
			rv := IDcount[rk]
			unionArr := GetUnionArr(v, rv)

			if len(unionArr) == 0 {
				nounionNum += 1
				continue
			}
			symIDC := unionArr[0]
			rcsymID := symIDC.ID + "_RC"
			symUnionArr := GetUnionArr(IDcount[symIDC.ID], IDcount[rcsymID])
			if len(symUnionArr) == 0 {
				nounionNum += 1
				continue
			}
			if symUnionArr[0].ID == k {
				symmetryMatch += 1
			}
			// sz := readLenMap[symIDC.ID]
			// segNum := (sz - SegmentLen) / StepLen
			fmt.Printf("%s\t%v\n", k, unionArr)
			fmt.Printf("%s\t%v\n\n", symIDC.ID, symUnionArr)
		}
	}

	fmt.Printf("symmetryMatch: %d, totalNum: %d, percent: %f\n", symmetryMatch, totalNum, float64(symmetryMatch)/float64(totalNum))

	return
}

func DCT(cmd cli.Command) {
	// C.jpgcvDCT()
	// C.ReadsDCT()
	// double* ReadsDCT(double data[], const int length) {
	faName := cmd.Parent().Flag("f").String()
	prefix := cmd.Parent().Flag("p").String()
	arg, err := strconv.Atoi(cmd.Flag("Depth").String())
	if err != nil {
		log.Fatalf("[DCT] Depth Get() error: %v\n", err)
	} else {
		Depth = arg
	}
	SVDAfn := prefix + ".svdA"
	SVDAfp, err := os.Create(SVDAfn)
	// defer SVDAfp.Close()
	if err != nil {
		log.Fatalf("[DCT] create %s err : %v\n", SVDAfn, err)
	}

	fafp, err := os.Open(faName)
	if err != nil {
		log.Fatal("[DCT]open fa input file failed")
	}
	defer fafp.Close()
	fagzip, err := gzip.NewReader(fafp)
	far := fasta.NewReader(fagzip, linear.NewSeq("", nil, alphabet.DNA))
	// var dctMatArr []float64
	var rows, cols int = 5, 5
	var topSigmaNum int = 25
	var segment_len int = SegmentLen
	var step_len int = StepLen
	var svdVectArr SvdVectArr
	var vectNum = 0
	var lastDctArr []float64
	readLenMap := make(map[string]int)
	// min := 100000

	for {
		if s, err := far.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("failed read %v:%s\n", faName, err)
			}
		} else {
			// fmt.Printf("ID: %v\n", s.CloneAnnotation().ID)
			t := s.(*linear.Seq)
			tl := t.Slice().(alphabet.Letters)
			// skip the length smaller than segment_len
			if len(tl) < segment_len {
				continue
			}
			digitArr := make([]DigitSrc, 2)
			digitArr[0].ID = s.CloneAnnotation().ID
			digitArr[0].src = DNA2Pixel(tl)
			digitArr[1].ID = digitArr[0].ID + "_RC"
			digitArr[1].src = DNA2PixelRev(tl)
			readLenMap[digitArr[0].ID] = len(tl)

			for _, dsrc := range digitArr {
				for i := 0; i < len(dsrc.src)-segment_len; i += step_len {
					var svdVect SvdVect
					svdVect.ID = dsrc.ID
					svdVect.Position = i
					svdVectArr = append(svdVectArr, svdVect)
					// if len(srcarr) < min {
					// 	min = len(srcarr)
					// }
					// srcarr = AdjustSrcArr(srcarr)
					segment_arr := dsrc.src[i : i+segment_len]
					fmt.Printf("array len: %d, segment_arr: %d:%f, %d:%f\n", len(segment_arr), 0, segment_arr[0], len(segment_arr)-1, segment_arr[len(segment_arr)-1])
					// csrcarr := C.NewdoubleArray(C.int(len(srcarr)))
					// defer C.free(unsafe.Pointer(csrcarr))
					// fmt.Printf("Transform before:%v\n", csrcarr)
					// Gofloat64Arr2CdoubleArr(srcarr, unsafe.Pointer(csrcarr))
					// fmt.Printf("Transform after:%v\n", csrcarr)
					dctRet := C.ReadsDCT(&segment_arr[0], C.int(len(segment_arr)))
					cdctarr := dctRet.cdctarr
					dst := dctRet.dst
					// defer C.free(unsafe.Pointer(cdctarr))
					dctarr := CdoubleArr2Gofloat64Arr(unsafe.Pointer(cdctarr), len(segment_arr))
					// fmt.Printf("dctarr: %v\n", cdctarr)
					// PrintCdctarr(cdctarr, len(srcarr))
					primedct := GetPrimeDCTArr(dctarr, len(segment_arr), 1, rows*cols)
					if len(primedct) != rows*cols {
						log.Fatalf("len(transdct):%d	!= rows*cols:%d\n", len(primedct), rows*cols)
					}
					// write prime dct to output file
					PrintPrimeDCT(svdVect, primedct, 1, rows*cols)
					if vectNum > 0 {
						CompareDctArr(lastDctArr, primedct)
					}
					for j, e := range primedct {
						// sl := fmt.Sprintf("%d\t%d\t%f\n", vectNum, j, e)
						// SVDAfp.WriteString(sl)
						fmt.Fprintf(SVDAfp, "%d\t%d\t%f\n", vectNum, j, e)
					}
					vectNum += 1

					// dctMatArr = append(dctMatArr, primedct...)
					// var a DCTArr
					// a.dct = dctarr
					// dctMat = append(dctMat, a)
					// C.free(unsafe.Pointer(cdctarr))
					C.cvReleaseMat(&dst)
					lastDctArr = primedct
				}

			}

		}
	}

	// compute SVD matrix
	SVDAfp.Close()
	// call system spark application
	sparkcmd := exec.Command("/db/software/spark-1.3.1-bin-hadoop2.6/bin/spark-submit", "--class", "SVD",
		"--driver-memory", "6G", "--executor-memory", "12G", "--driver-cores", "2",
		"/db/goworkspace/src/pkc/test/svd-scala-spark/target/scala-2.10/svd_2.10-1.0.jar", SVDAfn, strconv.Itoa(topSigmaNum))
	sparkstdout, _ := os.Create("spark.stdout")
	defer sparkstdout.Close()
	sparkcmd.Stdout = sparkstdout
	sparkstderr, _ := os.Create("spark.stderr")
	defer sparkstderr.Close()
	sparkcmd.Stderr = sparkstderr
	if err := sparkcmd.Run(); err != nil {
		log.Fatalf("[DCT]spark Command called error: %v\n", err)
	}

	// read U files
	Ufiles, err := CheckSuccAndListDir("U")
	if err != nil {
		log.Fatal(err)
	}
	for _, fn := range Ufiles {
		fp, err := os.Open(fn)
		if err != nil {
			log.Fatal(err)
		}
		defer fp.Close()
		buffp := bufio.NewReader(fp)
		eof := false
		for !eof {
			line, err := buffp.ReadString('\n')
			if err != nil {
				if err == io.EOF {
					err = nil
					eof = true
					continue
				} else {
					log.Fatal(err)
				}
			}

			// ts := strings.TrimSpace(line)
			fds := strings.Fields(line)
			if len(fds) != topSigmaNum+1 {
				log.Fatalf("len(fds):%d != %d\n", len(fds), topSigmaNum+1)
			}
			id, err := strconv.Atoi(fds[0])
			svdVectArr[id].LeftSV = StrArr2float64Arr(fds[1:])
		}
	}

	// Check and read Sigma file
	{
		Sigmafiles, err := CheckSuccAndListDir("Sigma")
		if err != nil {
			log.Fatal(err)
		}

		if len(Sigmafiles) != 1 {
			log.Fatalf("Sigmafile file must be one file\n")
		}

		sigmafp, err := os.Open(Sigmafiles[0])
		if err != nil {
			log.Fatal(err)
		}
		defer sigmafp.Close()
		buffp := bufio.NewReader(sigmafp)
		eof := false
		for !eof {
			line, err := buffp.ReadString('\n')
			if err != nil {
				if err == io.EOF {
					err = nil
					eof = true
					continue
				} else {
					log.Fatal(err)
				}
			}

			fds := strings.Fields(line)
			if len(fds) != 1 {
				log.Fatalf("len(fds):%d != %d\n", len(fds), 1)
			}
			fe, _ := strconv.ParseFloat(fds[0], 64)
			TopSigmaArr = append(TopSigmaArr, fe)
		}
	}

	// vectNum :=len(dctMatArr) / (rows * cols)
	// fmt.Printf("[SVD]len(dctMatArr): %d, vectNum: %v, len(svdVectArr): %d\n", len(dctMatArr), vectNum, len(svdVectArr))
	// // A := smatrix.MakeDenseMatrix(dctMatArr, vectNum, rows*cols)

	// A := mat64.NewDense(vectNum, rows*cols, dctMatArr)
	// fmt.Printf("[SVD]initial matrix A: %v\n", A)
	// epsilon := math.Pow(2, -52.0)
	// small := math.Pow(2, -966.0)
	// svdf := mat64.SVD(A, epsilon, small, true, false)

	// // U, Σ, V, err := A.SVD()
	// if err != nil {
	// 	log.Fatal(err)
	// }

	// // U, Σ, _, _ := A.SVD()
	// // svM := C.SVD(&dctMatArr[0], C.int(vectNum), C.int(rows*cols))
	// // C.free(unsafe.Pointer(svM.vt))
	// // fmt.Printf("\n[SVD]svM: %v, svM.s: %v, svM.u: %v\n", svM, svM.s, svM.u)
	// // os.Stdout.Sync()
	// // singluarV := CdoubleArr2Gofloat64Arr(unsafe.Pointer(svM.s), rows*cols)
	// // C.free(unsafe.Pointer(svM.s))
	// fmt.Printf("[SVD]Left singluar U: %v\n", svdf)
	// fmt.Printf("[SVD]singluar value Sigma: %v\n", svdf.Sigma)
	// // fmt.Printf("[SVD]Left singluar V: %v\n", V)
	// // os.Stdout.Sync()
	// // uarrays := U.Arrays()
	// for i := 0; i < vectNum; i++ {
	// 	svdVectArr[i].LeftSV = svdf.U.Row(nil, i)
	// 	// fmt.Printf("[SVD]LeftSV: %v\n", svdVectArr[i])
	// }

	// // sort left Singluar Vectors
	// fmt.Printf("After sort Left Singluar Vectors\n")
	// os.Stdout.Sync()
	sort.Sort(svdVectArr)
	for _, item := range svdVectArr {
		fmt.Printf("%v/%v: %10.2v\n", item.ID, item.Position, item.LeftSV)
	}

	// get every reads nearst same read order
	// NearestReadInfoMap := ConstructNearestReadRelation(svdVectArr, readLenMap)
	ConstructNearestReadRelation(svdVectArr, readLenMap)
}
