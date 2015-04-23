package discretecosinetransform

// #cgo CFLAGS: -I/usr/local/include
// #cgo LDFLAGS: -L ./ -lcv -lcxcore -ljpeg -lpng -lhighgui
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

float* ReadsDCT(float data[], const int length) {
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
	"github.com/jwaldrip/odin/cli"
)

func DCT(c cli.Command) {
	// C.jpgcvDCT()
	C.ReadsDCT()
	// float* ReadsDCT(float data[], const int length) {
}
