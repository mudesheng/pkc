package discretecosinetransform

// #cgo CFLAGS: -I/usr/local/include -std=c99
// #cgo LDFLAGS: -L ./ -lm -lcv -lcxcore -ljpeg -lpng -lhighgui
// #include <cv.h>
// #include <stdio.h>
// #include <math.h>
// #include <highgui.h>
// #include <cxcore.h>
/*
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

float* NewFloatArray(int length) {
  float *arr = calloc(length, sizeof(float));
  return arr;
}

inline void SetFloatArray(float *arr, int pos, float val) {
  arr[pos] = val;
}

float* ReadsDCT(float *data, const int length) {
  printf("data pointer address: %p\t", data);
  printf("length: %d, clang arr:\t", length);
  printf("%d:%f, %d:%f\n", 0, data[0], length-1, data[length-1]);
  // for( int i = 0; i < length; i++) {
  //   printf("%f ");
  // }
  printf("\n");
  double row = sqrt(length);
  if(row * row != length  && (int)row % 2 == 0) {
    fprintf(stderr, "[ReadsDCT]the input data matrix rows and columns must be equaled and even, exit...\n");
      exit(1);
  }

  CvMat src32 ;
  CvMat *dst32 = NULL;
  cvInitMatHeader(&src32, row,row ,CV_32FC1, data, CV_AUTOSTEP);
  cvDCT(&src32,dst32,CV_DXT_FORWARD);

  return (float*)dst32->data.fl;

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
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"unsafe"

	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
)

type DCTArr struct {
	dct []float32
}

var DNA2Num = [5]float32{0, 63, 127, 191, 255}

func DNA2Pixel(tl alphabet.Letters) (srcarr []float32) {
	for _, e := range tl {
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

	return
}

func AdjustSrcArr(dctarr []float32) []float32 {
	sq := math.Sqrt(float64(len(dctarr)))
	if sq-math.Floor(sq) != 0 || sq/2 != 0 {
		size := math.Floor(sq) + 1
		if size/2 != 0 {
			size += 1
		}
		size *= size
		for i := len(dctarr); i < int(size); i++ {
			dctarr = append(dctarr, DNA2Num[0])
		}
	}

	return dctarr
}

func GoFloat32Arr2CFloatArr(srcarr []float32, csrcarr unsafe.Pointer) {
	p := uintptr(csrcarr)
	for _, e := range srcarr {
		// *(*C.float)(unsafe.Pointer(&p)) = C.float(e)
		*(*C.float)(unsafe.Pointer(&p)) = C.float(e)
		fmt.Printf("%v\t", p)
		fmt.Printf("%f\t", *(*float32)(csrcarr))
		p += unsafe.Sizeof(e)
		fmt.Printf("%v\t", p)
	}
	fmt.Printf("%f\t", *(*float32)(csrcarr))
}

func CFloatArr2GoFloat32Arr(cdctarr unsafe.Pointer, size int) (dctarr []float32) {
	p := uintptr(cdctarr)

	for i := 0; i < size; i++ {
		j := *(*float32)(unsafe.Pointer(p))
		dctarr = append(dctarr, j)
		p += unsafe.Sizeof(j)
	}

	return dctarr
}

func DCT(cmd cli.Command) {
	// C.jpgcvDCT()
	// C.ReadsDCT()
	// float* ReadsDCT(float data[], const int length) {
	faName := cmd.Parent().Flag("f").String()
	fafp, err := os.Open(faName)
	if err != nil {
		log.Fatal("[DCT]open fa input file failed")
	}
	defer fafp.Close()
	fagzip, err := gzip.NewReader(fafp)
	far := fasta.NewReader(fagzip, linear.NewSeq("", nil, alphabet.DNA))
	var dctMat []DCTArr

	for {
		if s, err := far.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("failed read %v:%s\n", faName, err)
			}
		} else {
			t := s.(*linear.Seq)
			tl := t.Slice().(alphabet.Letters)
			srcarr := DNA2Pixel(tl)
			AdjustSrcArr(srcarr)
			fmt.Printf("array len: %d, srcarr: %d:%f, %d:%f\n", len(srcarr), 0, srcarr[0], len(srcarr)-1, srcarr[len(srcarr)-1])
			csrcarr := C.NewFloatArray(C.int(len(srcarr)))
			defer C.free(unsafe.Pointer(csrcarr))
			fmt.Printf("Transform before:%v\n", csrcarr)
			GoFloat32Arr2CFloatArr(srcarr, unsafe.Pointer(csrcarr))
			fmt.Printf("Transform after:%v\n", csrcarr)
			cdctarr := C.ReadsDCT(csrcarr, C.int(len(srcarr)))
			dctarr := CFloatArr2GoFloat32Arr(unsafe.Pointer(cdctarr), len(srcarr))
			fmt.Printf("dctarr: %v\n", dctarr)
			var a DCTArr
			a.dct = dctarr
			dctMat = append(dctMat, a)

		}
	}
}
