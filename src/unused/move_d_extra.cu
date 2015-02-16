#if UN_COMPILE
  // perturbation depends on the step number
  if (step%2 == 0){
    if (bidx < 6) {
      movematrix_new[bidx] = 1.0f;
    }
    if (bidx > 2 && bidx < 6) {    
      movematrix_new[bidx] = -1.0f;
    }
    if (bidx < 6) {
      movematrix_new[bidx] += mylig->movematrix_old[bidx];
      mylig->movematrix_new[bidx] = movematrix_new[bidx];
    }
  }
  else
    {
      if (bidx < 6) {
	movematrix_new[bidx] = -1.0f;
      }
      if (bidx > 2 && bidx < 6) {    
	movematrix_new[bidx] = 1.0f;
      }
      if (bidx < 6) {
	movematrix_new[bidx] += mylig->movematrix_old[bidx];
	mylig->movematrix_new[bidx] = movematrix_new[bidx];
      }
    }
#endif
 

#if UN_COMPILE
  // only translational vector applied
  if (bidx < 6) {
    movematrix_new[bidx] = 1.0f;
  }
  if (bidx > 2 && bidx < 6) {
    movematrix_new[bidx] = 0.0f;
  }
 if (bidx < 6) {
    movematrix_new[bidx] += mylig->movematrix_old[bidx];
    mylig->movematrix_new[bidx] = movematrix_new[bidx];
  }
#endif
 
#if UN_COMPILE
  // only rotational vector applied
  if (bidx < 6) {
    movematrix_new[bidx] = 1.0f;
  }
  if (bidx < 3) {
    movematrix_new[bidx] = 0.0f;
  }
 if (bidx < 6) {
    movematrix_new[bidx] += mylig->movematrix_old[bidx];
    mylig->movematrix_new[bidx] = movematrix_new[bidx];
  }
#endif



