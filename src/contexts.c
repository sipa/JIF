#include "contexts.h"

void init_context_leaf_counts(context_tree *first) {
  first->property = 0;
  first->data.context->count = 0;
  first->data.context->bits = 0;
  for (int i=0 ; i < MAX_NB_PROPERTIES ; i++) {
          first->data.context->sum[i] = 0;
          first->data.context->splitbits[i] = 1;
  }
}
void init_context_leaf_splitcontexts(context_tree *first) {
  init_context_leaf_counts(first);
  for (int i=0 ; i < MAX_NB_PROPERTIES*2 ; i++) {
     symb_chs_init( &(first->data.context->splitcontext[i]));
  }
}
void init_context_leaf(context_tree *first) {
  first->data.context = malloc(sizeof(context_leaf));
  symb_chs_init( &(first->data.context->context));
  init_context_leaf_splitcontexts(first);
}

void init_context_leaf_counts_bit(context_tree_bit *first) {
  first->property = 0;
  first->data.context->count = 0;
  first->data.context->bits = 0;
  for (int i=0 ; i < NB_PROPERTIES_S ; i++) {
          first->data.context->sum[i] = 0;
          first->data.context->splitbits[i] = 1;
  }
}
void init_context_leaf_splitcontexts_bit(context_tree_bit *first) {
  init_context_leaf_counts_bit(first);
  for (int i=0 ; i < NB_PROPERTIES_S*2 ; i++) {
     chs_init( &(first->data.context->splitcontext[i]),2048);
  }
}
void init_context_leaf_bit(context_tree_bit *first) {
  first->data.context = malloc(sizeof(context_leaf_bit));
  chs_init( &(first->data.context->context), 2048);
  init_context_leaf_splitcontexts_bit(first);
}


void ctx_init(ctx_t *ctx, int cutoff, int phase) {

  // initialize split chances
  ctx->bitsSplit[3][0]=0;
  ctx->symbSplit[3][0]=0;
  ctx->bitsSplit[3][1]=0;
  ctx->symbSplit[3][1]=0;

  ctx->bitsContextTreeData[0]=0;   
  ctx->bitsContextTreeData[1]=0;   
  ctx->bitsContextTreeData[2]=0;   

//  ctx->bitsColorData=0;
  ctx->bitsColorData[0] = 0;
  ctx->bitsColorData[1] = 0;
  ctx->bitsColorData[2] = 0;
  ctx->bitsColorData_Y=0;
  ctx->bitsColorData_count=0;



  for (int c = 0; c < 3; c++) {
    ctx->bitsSplit[c][0]=0;
    ctx->symbSplit[c][0]=0;
    ctx->bitsSplit[c][1]=0;
    ctx->symbSplit[c][1]=0;
    ctx->bitsInterpol[c]=0;
    ctx->symbInterpol[c]=0;
    ctx->bitsPixeldata[c][0]=0;
    ctx->symbPixeldata[c][0]=0;
    ctx->bitsPixeldata[c][1]=0;
    ctx->symbPixeldata[c][1]=0;
/*
    for (int h = 0; h < QZ_SIZE_S; h++) 
    for (int i = 0; i < QZ_DEPTH_S; i++) 
    for (int j = 0; j < 4; j++) 
    for (int n = 0; n < QZ_BVAR_S; n++) {
         chs_init(&ctx->splitContA[h][i][c][j][n],clamp(50*(c+5)*(QZ_SIZE_S-h+QZ_DEPTH_S-i)/(QZ_SIZE_S+QZ_DEPTH_S) + 25*(c) + 50*(c+5)*(QZ_BVAR_S-n)/QZ_BVAR_S,cutoff,4096-cutoff));
//         chs_init(&ctx->splitContA[h][i][c][j][n],2048);
        }

    for (int h = 0; h < QZ_SIZE_SL; h++) 
    for (int i = 0; i < QZ_DEPTH_SL; i++) 
    for (int j = 0; j < 7; j++) 
    for (int n = 0; n < QZ_BVAR_S; n++) {
         for (int m = 0; m < QZ_LRDIFF_S; m++) {
            chs_init(&ctx->splitContL[h][i][c][j][n][m],clamp(50*(c+5)*(QZ_SIZE_SL-h+QZ_DEPTH_SL-i)/(QZ_SIZE_SL+QZ_DEPTH_SL) + 25*(c) + 50*(c+5)*(QZ_BVAR_S-n)/QZ_BVAR_S,cutoff,4096-cutoff));
            //2048);
         }
    }
*/
    for (int i = 0; i < QZ_DEPTH_S; i++) 
    for (int j = 0; j < MAX_INTERPOLATIONS; j++) 
    for (int k = 0; k < MAX_INTERPOLATIONS; k++) 
          chs_init(&ctx->interpolationMethod[i][c][j][k],2048);
  }

  context_tree *firstY; 
  context_tree *firstQ; 
  context_tree *firstI; 

  context_tree_bit *firstSY;
  context_tree_bit *firstSQ;
  context_tree_bit *firstSI;

  if (phase != 2) {
  ctx->treesize[0][0] = 0;
  ctx->treesize[0][1] = 0;
  ctx->treesize[0][2] = 0;
  firstY = &ctx->diff[0][0];
  firstQ = &ctx->diff[1][0];
  firstI = &ctx->diff[2][0];
  init_context_leaf(firstY);
  init_context_leaf(firstI);
  init_context_leaf(firstQ);
  }

#if FULL_BORDERS
  ctx->treesize[3][0] = 0;
  ctx->treesize[3][1] = 0;
  ctx->treesize[3][2] = 0;
  firstY = &ctx->diff1[0][0];
  firstQ = &ctx->diff1[1][0];
  firstI = &ctx->diff1[2][0];
  init_context_leaf(firstY);
  init_context_leaf(firstI);
  init_context_leaf(firstQ);

  ctx->treesize[1][0] = 0;
  ctx->treesize[1][1] = 0;
  ctx->treesize[1][2] = 0;
  firstSY = &ctx->splitContL[0][0];
  firstSQ = &ctx->splitContL[1][0];
  firstSI = &ctx->splitContL[2][0];
  init_context_leaf_bit(firstSY);
  init_context_leaf_bit(firstSI);
  init_context_leaf_bit(firstSQ);

#endif


  ctx->treesize[2][0] = 0;
  ctx->treesize[2][1] = 0;
  ctx->treesize[2][2] = 0;
  firstSY = &ctx->splitContA[0][0];
  firstSQ = &ctx->splitContA[1][0];
  firstSI = &ctx->splitContA[2][0];
  init_context_leaf_bit(firstSY);
  init_context_leaf_bit(firstSI);
  init_context_leaf_bit(firstSQ);


  symb_chs_init(&ctx->methods);

  for (int c = 0; c < 3; c++) {
    symb_chs_init(&ctx->r_min[c]);
    symb_chs_init(&ctx->r_max[c]);
    symb_chs_init(&ctx->colors[c]);
  }
  symb_chs_init(&ctx->colors_Y);
  symb_chs_init(&ctx->colors_count0);
  symb_chs_init(&ctx->colors_count);
  symb_chs_init(&ctx->context_tree_property);
  symb_chs_init(&ctx->context_tree_count1);
  symb_chs_init(&ctx->context_tree_count2);
}

static const char *context_properties[] = {
  "I-Diff",
  "Y-Diff",
  "Y-Guess",
  "I-Guess",
  "Q-Guess",
  "Size",
  "LR-diff",
  "Variance",
  "Border variance",
  "Border diff expected",
  "Split-direction",
// not a very good idea to add these:
//  "X",
//  "Y"
  NULL
};
#ifdef DEBUGMODE
static const char *context_properties_1[] = {
  "I-Diff",
  "Y-Diff",
  "Y-Guess",
  "I-Guess",
  "Q-Guess",
  "1-5 Diff",
  "3-7 Diff",
  "2-6 Diff",
  "8-4 Diff",
  "Variance",
  "Border variance",
  NULL
};
#if SPLIT_CHANNELS_SEPARATELY     
static const char *context_properties_bit[] = {
  "Size",
  "LR-diff / Border diff",
  "Border variance",
  "Other splits",
  "Depth-Avg",
  NULL
};
#else
static const char *context_properties_bit[] = {
  "Size",
  "Border var Y",
  "Border var I+Q",
  "Other splits",
  "Depth-Avg",
  NULL
};
#endif
#endif


void output_context_subtree(context_tree *tree, coder_t *coder, int treenumber, int c, int pos ) {
        int p=tree[pos].property;
        int nb_properties;
        if (coder->traversal_method == 0) {
                nb_properties = (treenumber == 0 ? NB_PROPERTIES : NB_PROPERTIES_1);
        } else {
                nb_properties = NB_PROPERTIES_FFV1;
        }

//        fprintf(stderr,"Entering pos %i, property %i\n",pos,p);

        symb_put_int(&coder->coder,&coder->ctx->context_tree_property,p,0,nb_properties,&(coder->ctx->bitsContextTreeData[c]));
        
        
        if (p == 0) return;

//        int splitval = tree[pos].data.node.splitval;
                //average(tree[pos].data.node.sum, tree[pos].data.node.count);

//        symb_put_simple_int(&coder->coder,splitval,lower[p-1],upper[p-1],&(coder->ctx->bitsContextTreeData[c]));

        int count = tree[pos].data.node.count;
        if (count > CONTEXT_SPLIT_MAX_SIZEP) {
                fprintf(stderr,"Warning: context node count %i > %i, clamping value.\n",count, CONTEXT_SPLIT_MAX_SIZE);
                count = CONTEXT_SPLIT_MAX_SIZEP;
        }
        symb_put_big_int(&coder->coder,&coder->ctx->context_tree_count1,&coder->ctx->context_tree_count2,count,CONTEXT_SPLIT_MIN_SIZEP,CONTEXT_SPLIT_MAX_SIZEP,&(coder->ctx->bitsContextTreeData[c]));
        

//        fprintf(stderr,"Pos: %i, property %i (%s), count %i\n", pos, tree[pos].property-1, context_properties[tree[pos].property-1],count);
//        fprintf(stderr,"Pos: %i, property %i (%s), splitval %i, count %i\n", pos, tree[pos].property-1, context_properties[tree[pos].property-1], splitval, count);
        
        int pos1 = tree[pos].data.node.branch;
        output_context_subtree(tree, coder, treenumber, c, pos1);

        int pos2 = tree[pos].data.node.branch+1;
        output_context_subtree(tree, coder, treenumber, c, pos2);
}

void output_context_tree(context_tree *tree, coder_t *coder, int treenumber, int c) {
        int nb_properties;
        if (coder->traversal_method == 0) {
                nb_properties = (treenumber == 0 ? NB_PROPERTIES : NB_PROPERTIES_1);
        } else {
                nb_properties = NB_PROPERTIES_FFV1;
        }
        output_context_subtree(tree, coder, treenumber, c, 0);
}


void output_context_subtree_p2(context_tree_p2 *tree, coder_t *coder, int treenumber, int c, int pos) {
        int nb_properties;
        if (coder->traversal_method == 0) {
                nb_properties = (treenumber == 0 ? NB_PROPERTIES : NB_PROPERTIES_1);
        } else {
                nb_properties = NB_PROPERTIES_FFV1;
        }
        int p=tree[pos].property;

//        fprintf(stderr,"Entering pos %i, property %i\n",pos,p);

        symb_put_int(&coder->coder,&coder->ctx->context_tree_property,p,0,nb_properties,&(coder->ctx->bitsContextTreeData[c]));
        
        
        if (p == 0) return;

//        int splitval = tree[pos].splitval;
                //average(tree[pos].data.node.sum, tree[pos].data.node.count);

//        symb_put_simple_int(&coder->coder,splitval,lower[p-1],upper[p-1],&(coder->ctx->bitsContextTreeData[c]));

        int count = tree[pos].count;
        symb_put_big_int(&coder->coder,&coder->ctx->context_tree_count1,&coder->ctx->context_tree_count2,count,CONTEXT_SPLIT_MIN_SIZEP,CONTEXT_SPLIT_MAX_SIZEP,&(coder->ctx->bitsContextTreeData[c]));
        

//        fprintf(stderr,"Pos: %i, property %i (%s), count %i\n", pos, tree[pos].property-1, context_properties[tree[pos].property-1], count);
        
        int pos1 = tree[pos].branch;
        output_context_subtree_p2(tree, coder, treenumber, c, pos1);

        int pos2 = tree[pos].branch+1;
        output_context_subtree_p2(tree, coder, treenumber, c, pos2);
}

void output_context_tree_p2(context_tree_p2 *tree, coder_t *coder, int treenumber, int c) {
        int nb_properties;
        if (coder->traversal_method == 0) {
                nb_properties = (treenumber == 0 ? NB_PROPERTIES : NB_PROPERTIES_1);
        } else {
                nb_properties = NB_PROPERTIES_FFV1;
        }
        output_context_subtree_p2(tree, coder, treenumber, c, 0);
}



void input_context_subtree_p2(context_tree_p2 *tree, coder_t *coder, int treenumber, int c, int pos, int *counter) {
        int nb_properties;
        if (coder->traversal_method == 0) {
                nb_properties = (treenumber == 0 ? NB_PROPERTIES : NB_PROPERTIES_1);
        } else {
                nb_properties = NB_PROPERTIES_FFV1;
        }
        int p = symb_get_int(&coder->coder,&coder->ctx->context_tree_property,0,nb_properties);
        tree[pos].property = p;
//        fprintf(stderr,"Entering pos %i, property %i\n",pos,p);
        if (p == 0) {
                return;
        }

//        int splitval = symb_get_simple_int(&coder->coder,lower[p-1],upper[p-1]);
//        tree[pos].splitval = splitval;
//        tree[pos].data.node.sum = splitval;
//        tree[pos].data.node.count = 1;
        tree[pos].branch = *counter;
        *counter = *counter + 2;

        int count = symb_get_big_int(&coder->coder,&coder->ctx->context_tree_count1,&coder->ctx->context_tree_count2,CONTEXT_SPLIT_MIN_SIZEP,CONTEXT_SPLIT_MAX_SIZEP);
        tree[pos].count = count;
        tree[pos].nb = 0;
        tree[pos].sum = 0;
        tree[pos].splitval = 0;
        
//        fprintf(stderr,"Pos: %i, property %i (%s), count %i\n", pos, tree[pos].property-1, context_properties[tree[pos].property-1],count);
        
        int pos1 = tree[pos].branch;
        input_context_subtree_p2(tree, coder, treenumber, c, pos1, counter);

        int pos2 = tree[pos].branch+1;
        input_context_subtree_p2(tree, coder, treenumber, c, pos2, counter);
}

void input_context_tree_p2(context_tree_p2 *tree, coder_t *coder, int treenumber, int c) {
        int nb_properties;
        if (coder->traversal_method == 0) {
                nb_properties = (treenumber == 0 ? NB_PROPERTIES : NB_PROPERTIES_1);
        } else {
                nb_properties = NB_PROPERTIES_FFV1;
        }
        int counter = 1;
        input_context_subtree_p2(tree, coder, treenumber, c, 0, &counter);

        tree[0].context = malloc(sizeof(symb_chs_t));
        symb_chs_init(tree[0].context);

        // set treesize to max to avoid further tree construction
//        int max_contexts = (treenumber == 0 ? MAX_CONTEXTS : MAX_CONTEXTS_1);
//        coder->ctx->treesize[treenumber][c] = max_contexts;
//        coder->ctx->treesize[treenumber][c] = counter;
}



//const int upper_init[NB_PROPERTIES] = {QZ_ID-1,QZ_YD-1,QZ_YG-1,QZ_IG-1,QZ_QG-1,QZ_SIZE-1,QZ_LRDIFF-1,QZ_VAR-1,QZ_BVAR-1};

inline static symb_chs_t* find_context_p2(context_tree_p2 *tree,int *properties,coder_t *coder,int treenumber, int c) {
        int pos=0;
        int copy_parent = -1;
        while(1) {
                int p=tree[pos].property;
                if (p == 0) break;
                // int splitval = average(tree[pos].data.node.sum, tree[pos].data.node.count);

                if (tree[pos].count > 0) {
                    tree[pos].nb++;
                    tree[pos].sum += properties[p-1];
                    tree[pos].count--;
                    if (tree[pos].count > 0) { 
                       break; 
                    } else { 
                       copy_parent = pos; 
                       tree[copy_parent].splitval = average(tree[copy_parent].sum,tree[copy_parent].nb);
                    }
                }
                int splitval = tree[pos].splitval;
                if (properties[p-1] > splitval) {
                        pos = tree[pos].branch;
                } else {
                        pos = tree[pos].branch+1;
                }
                if (copy_parent >= 0) break;
        }

          if (copy_parent >= 0) {
               int child = tree[copy_parent].branch;
//               symb_cp(tree[copy_parent].context,tree[child].context);
               tree[child].context = tree[copy_parent].context;
               tree[child+1].context = malloc(sizeof(symb_chs_t));
               symb_cp(tree[copy_parent].context,tree[child+1].context);
               tree[pos].sum += properties[tree[pos].property-1];
               tree[pos].count--;
               tree[pos].nb++;
//               fprintf(stderr,"2Channel %i: node %i, property %i, splitval %i = %lli / %i\n",c,copy_parent, tree[copy_parent].property, tree[copy_parent].splitval, tree[copy_parent].sum, tree[copy_parent].nb);
          }          
//          fprintf(stderr,"Channel %i, pos %i (remaining %i)\n",c,pos,tree[pos].count);
          return tree[pos].context;
        
}




inline context_leaf* find_context(context_tree *tree,int *properties,coder_t *coder,int treenumber, int c, int *vcontext) {
        int pos=0;
        int nb_properties;
        if (coder->traversal_method == 0) {
#if FULL_BORDERS
                nb_properties = (treenumber == 0 ? NB_PROPERTIES : NB_PROPERTIES_1);
#else
                nb_properties = NB_PROPERTIES; 
#endif
        } else {
                nb_properties = NB_PROPERTIES_FFV1;
        }
#if FULL_BORDERS
        int max_contexts = (treenumber == 0 ? MAX_CONTEXTS : MAX_CONTEXTS_1);
#else
        const int max_contexts = MAX_CONTEXTS;
#endif
        int upper[MAX_NB_PROPERTIES];
        int lower[MAX_NB_PROPERTIES] = {};

//        for (int i = 0; i<nb_properties; i++) lower[i]=0;
        
        upper[0]=QZ_ID-1;
        upper[1]=QZ_YD-1;
        upper[2]=QZ_YG-1;
        upper[3]=QZ_IG-1;
        upper[4]=QZ_QG-1;

        if (coder->traversal_method == 0) {     // splitting
#if FULL_BORDERS
          if (treenumber == 0) {
#endif
                upper[5]=QZ_SIZE-1;
//                upper[6]=DEPTH_LEVELS-1;
                upper[6]=QZ_LRDIFF-1;
                upper[7]=QZ_VAR-1;
                upper[8]=QZ_BVAR-1;
                upper[9]=QZ_BDIFFE-1;
                upper[10]=1;
//                upper[9]=MAX_X;
//                upper[10]=MAX_Y;
#if FULL_BORDERS
          } else {
                upper[5]=QZ_LRDIFF-1;
                upper[6]=QZ_LRDIFF-1;
                upper[7]=QZ_LRDIFF-1;
                upper[8]=QZ_LRDIFF-1;
                upper[9]=QZ_VAR-1;
                upper[10]=QZ_BVAR-1;
          }                
#endif
        } else {                                // FFV1
          for (int i = 5; i<nb_properties; i++) upper[i]=QZ_LRDIFF-1;
        }

        int max_size = CONTEXT_SPLIT_MIN_SIZE;

        while(1) {
                int p=tree[pos].property;
                if (p == 0) break;
                // int splitval = average(tree[pos].data.node.sum, tree[pos].data.node.count);
                int splitval = tree[pos].data.node.splitval;

// uncomment for adaptive splitvals (seems to be bad idea)
//                tree[pos].data.node.sum += properties[p-1];
//                tree[pos].data.node.count++;

                if (properties[p-1] > splitval) {
                        lower[p-1] = splitval+1;
                        pos = tree[pos].data.node.branch;
                } else {
                        upper[p-1] = splitval;
                        pos = tree[pos].data.node.branch+1;
                }
#if MAX_CONTEXT_SIZE_DEPTH_FACTOR > 1
                max_size = MAX_CONTEXT_SIZE_DEPTH_FACTOR * max_size;
#endif                
        }
        context_leaf *current = tree[pos].data.context;
        
        int return_pos = pos;
        current->count++;
        uint64_t min_bits=0;
        int best_split_property=0;
//        fprintf(stderr,"Current context: %i \n",pos);
        for (int i=2-c ; i<nb_properties ; i++) {
              uint64_t bits = current->splitbits[i];
//                int range = upper[i]-lower[i];
#if (IGNORE_SPLIT_OPTIONS==1)
                if (bits>0 && range>0 && (current->count < CONTEXT_SPLIT_IGNORE_COUNT || bits < CONTEXT_SPLIT_IGNORE_FACTOR * current->bits)) {
#else
              if (bits>0) {
#endif
                if (upper[i]-lower[i] > 0) {
                  current->sum[i] += properties[i];
                  if (!min_bits || (bits < min_bits)) {
                        min_bits = bits;
                        best_split_property=i;
                  }
                } else {
                  current->splitbits[i]=0;
//                  fprintf(stderr,"Empty range (%i..%i) for property %i (%s)\n",lower[i],upper[i],i,context_properties[i]);
                }
              }
        }
        if (current->count == CONTEXT_VIRTUAL_OUTPUT_DELAY)  {
             for (int i=2-c ; i<nb_properties ; i++) {
                symb_cp(&tree[pos].data.context->context,&tree[pos].data.context->splitcontext[i]);
                symb_cp(&tree[pos].data.context->context,&tree[pos].data.context->splitcontext[i+nb_properties]);
             }
        }
        
        if ( //current->count >= CONTEXT_SPLIT_MAX_SIZE ||
(                current->count > CONTEXT_VIRTUAL_OUTPUT_DELAY && current->count > max_size 
                && coder->ctx->treesize[treenumber][c] < max_contexts 
                && current->bits > CONTEXT_SPLIT_IMPROVE_CONSTANT + CONTEXT_SPLIT_IMPROVE_FACTOR * current->splitbits[best_split_property])) {
//                if (treenumber == 3)                fprintf(stderr,"Splitting on property %s\n", context_properties_1[best_split_property]);
//                fprintf(stderr,"Splitting on property %s, splitval %i (range:%i..%i), count:%u , bits:%llu > %llu\n", context_properties_1[best_split_property], (int) (current->sum[best_split_property]/current->count),lower[best_split_property],upper[best_split_property], current->count, current->bits, current->splitbits[best_split_property]);
//                fprintf(stderr,"Splitting context %i on property %s, splitval %i (range:%i..%i), count:%u , bits:%llu > %llu\n", pos, context_properties[best_split_property], (int) (current->sum[best_split_property]/current->count),lower[best_split_property],upper[best_split_property], current->count, current->bits, current->splitbits[best_split_property]);
//                fprintf(stderr,"Splitting on property %s\n", context_properties[best_split_property]);
                int child1 = ++coder->ctx->treesize[treenumber][c];
                int child2 = ++coder->ctx->treesize[treenumber][c];
                if (properties[best_split_property] > average(current->sum[best_split_property],current->count)) {
                    return_pos = child1;
                } else {
                    return_pos = child2;
                }
                tree[child1].data.context = tree[pos].data.context;
                init_context_leaf(&tree[child2]);
                // set children contexts to virtual contexts 
                symb_cp(&tree[pos].data.context->splitcontext[nb_properties+best_split_property],&tree[child1].data.context->context);
                symb_cp(&tree[pos].data.context->splitcontext[best_split_property],&tree[child2].data.context->context);
/*
                for (int i=2-c ; i<nb_properties ; i++) {
                   if (i != best_split_property) {
                        symb_cp(&tree[pos].data.context->splitcontext[i],&tree[child1].data.context->splitcontext[i]);
                        symb_cp(&tree[pos].data.context->splitcontext[i],&tree[child2].data.context->splitcontext[i]);
                        symb_cp(&tree[pos].data.context->splitcontext[i+nb_properties],&tree[child1].data.context->splitcontext[i+nb_properties]);
                        symb_cp(&tree[pos].data.context->splitcontext[i+nb_properties],&tree[child2].data.context->splitcontext[i+nb_properties]);
                   } else {
                        symb_cp(&tree[pos].data.context->splitcontext[i+nb_properties],&tree[child1].data.context->splitcontext[i]);
                        symb_cp(&tree[pos].data.context->splitcontext[i+nb_properties],&tree[child1].data.context->splitcontext[i+nb_properties]);
                        symb_cp(&tree[pos].data.context->splitcontext[i],&tree[child2].data.context->splitcontext[i]);
                        symb_cp(&tree[pos].data.context->splitcontext[i],&tree[child2].data.context->splitcontext[i+nb_properties]);
                   }
                }
                */
                tree[pos].property = best_split_property+1;
//                tree[pos].data.node.sum = current->sum[best_split_property];
                tree[pos].data.node.count = current->count;
                tree[pos].data.node.splitval = average(current->sum[best_split_property],current->count);
                tree[pos].data.node.branch = child1;
 //               fprintf(stderr,"1Channel %i: node %i, property %i, splitval %i = %lli / %i\n",c,pos, tree[pos].property, tree[pos].data.node.splitval, current->sum[best_split_property],current->count);
                init_context_leaf_counts(&tree[child1]);
                tree[return_pos].data.context->count = 1;
                for (int i=2-c ; i<nb_properties ; i++) {
                  tree[return_pos].data.context->sum[i] += properties[i];
                }
        } else {
                if (current->bits >  CONTEXT_USE_IMPROVE_CONSTANT + CONTEXT_USE_IMPROVE_FACTOR * current->splitbits[best_split_property]) {
                        *vcontext = best_split_property;
                }
        }
        return tree[return_pos].data.context;
       
        
}





context_leaf_bit* find_context_bit(context_tree_bit *tree,int properties[NB_PROPERTIES_S],coder_t *coder,int treenumber, int c, int *vcontext) {
        int pos=0;
#if SPLIT_CHANNELS_SEPARATELY     
        int upper[NB_PROPERTIES_S]={QZ_SIZE-1,QZ_LRDIFF_S-1,QZ_BVAR_S-1,3}; //,DEPTH_LEVELS};
#else
        int upper[NB_PROPERTIES_S]={QZ_SIZE-1,QZ_BVAR_S-1,QZ_BVAR_S-1,3}; //,DEPTH_LEVELS};
#endif
        int lower[NB_PROPERTIES_S]={0,0,0,0};
        int max_size = CONTEXT_SPLIT_MIN_SIZE_S;

        while(1) {
                int p=tree[pos].property;
                if (p == 0) break;
//                int splitval = average(tree[pos].data.node.sum, tree[pos].data.node.count);
                int splitval = tree[pos].data.node.splitval;
// uncomment for adaptive splitvals (seems to be bad idea)
//                tree[pos].data.node.sum += properties[p-1];
//                tree[pos].data.node.count++;
                if (properties[p-1] > splitval) {
                        lower[p-1] = splitval+1;
                        pos = tree[pos].data.node.branch;
                } else {
                        upper[p-1] = splitval;
                        pos = tree[pos].data.node.branch+1;
                }
                max_size = MAX_CONTEXT_SIZE_DEPTH_FACTOR * max_size;
        }
        context_leaf_bit *current = tree[pos].data.context;
        int return_pos = pos;
        current->count++;
        uint64_t min_bits=0;
        int best_split_property=0;
//        fprintf(stderr,"Current context: %i \n",pos);
        for (int i=0 ; i<NB_PROPERTIES_S ; i++) {
                uint64_t bits = current->splitbits[i];
                int range = upper[i]-lower[i];
#if (IGNORE_SPLIT_OPTIONS==1)
                if (bits>0 && range>0 && (current->count < CONTEXT_SPLIT_IGNORE_COUNT || bits < CONTEXT_SPLIT_IGNORE_FACTOR * current->bits)) {
#else
                if (bits>0 && range>0 ) {
#endif
                  current->sum[i] += properties[i];
                } else {
                  current->splitbits[i]=0;
                }
                if (!min_bits || (range>0 && bits>0 && bits < min_bits)) {
                        min_bits = bits;
                        best_split_property=i;
                }
        }
        if (current->count > max_size && coder->ctx->treesize[treenumber][c] < MAX_CONTEXTS_S && current->bits > CONTEXT_SPLIT_IMPROVE_CONSTANT_S + CONTEXT_SPLIT_IMPROVE_FACTOR * current->splitbits[best_split_property]) {
//                if (best_split_property == 4) 
//                fprintf(stderr,"Splitting on property %s, splitval %i (range:%i..%i), count:%u , bits:%llu > %llu\n", context_properties_bit[best_split_property], (int) (current->sum[best_split_property]/current->count),lower[best_split_property],upper[best_split_property], current->count, current->bits, current->splitbits[best_split_property]);
//                fprintf(stderr,"Splitting context %i on property %s, splitval %i (range:%i..%i), count:%u , bits:%llu > %llu\n", pos, context_properties_bit[best_split_property], (int) (current->sum[best_split_property]/current->count),lower[best_split_property],upper[best_split_property], current->count, current->bits, current->splitbits[best_split_property]);
//                fprintf(stderr,"Splitting on property %s\n", context_properties_bit[best_split_property]);
                int child1 = ++coder->ctx->treesize[treenumber][c];
                int child2 = ++coder->ctx->treesize[treenumber][c];
                if (properties[best_split_property] > average(current->sum[best_split_property],current->count)) {
                    return_pos = child1;
                } else {
                    return_pos = child2;
                }
                tree[child1].data.context = tree[pos].data.context;
                init_context_leaf_bit(&tree[child2]);
                // set children contexts to virtual contexts 
                chs_cp(&tree[pos].data.context->splitcontext[NB_PROPERTIES_S+best_split_property],&tree[child1].data.context->context);
                chs_cp(&tree[pos].data.context->splitcontext[best_split_property],&tree[child2].data.context->context);
                for (int i=0 ; i<NB_PROPERTIES_S ; i++) {
                   if (i != best_split_property) {
                        chs_cp(&tree[pos].data.context->splitcontext[i],&tree[child1].data.context->splitcontext[i]);
                        chs_cp(&tree[pos].data.context->splitcontext[i],&tree[child2].data.context->splitcontext[i]);
                        chs_cp(&tree[pos].data.context->splitcontext[i+NB_PROPERTIES_S],&tree[child1].data.context->splitcontext[i+NB_PROPERTIES_S]);
                        chs_cp(&tree[pos].data.context->splitcontext[i+NB_PROPERTIES_S],&tree[child2].data.context->splitcontext[i+NB_PROPERTIES_S]);
                   } else {
                        chs_cp(&tree[pos].data.context->splitcontext[i+NB_PROPERTIES_S],&tree[child1].data.context->splitcontext[i]);
                        chs_cp(&tree[pos].data.context->splitcontext[i+NB_PROPERTIES_S],&tree[child1].data.context->splitcontext[i+NB_PROPERTIES_S]);
                        chs_cp(&tree[pos].data.context->splitcontext[i],&tree[child2].data.context->splitcontext[i]);
                        chs_cp(&tree[pos].data.context->splitcontext[i],&tree[child2].data.context->splitcontext[i+NB_PROPERTIES_S]);
                   }
                }
                tree[pos].property = best_split_property+1;
                tree[pos].data.node.splitval = average(current->sum[best_split_property],current->count);
//                tree[pos].data.node.sum = current->sum[best_split_property];
//                tree[pos].data.node.count = current->count;
                tree[pos].data.node.branch = child1;
                init_context_leaf_counts_bit(&tree[child1]);
                tree[return_pos].data.context->count = 1;
                for (int i=0 ; i<NB_PROPERTIES_S ; i++) {
                  tree[return_pos].data.context->sum[i] += properties[i];
                }
        } else {
                if (current->bits >  CONTEXT_USE_IMPROVE_CONSTANT_S + CONTEXT_USE_IMPROVE_FACTOR *current->splitbits[best_split_property]) {
                        *vcontext = best_split_property;
                }
        }

        return tree[return_pos].data.context;
}
void output_interpolation_methods(coder_t *encode) {
        symb_put_int(&encode->coder,
                     &encode->ctx->methods,
                     encode->nb_interpolation_methods-1,
                     0,
                     MAX_INTERPOLATIONS_EVER,
                     NULL
                    );
        for(int i=0; i<encode->nb_interpolation_methods; i++) {
                symb_put_int(&encode->coder,
                     &encode->ctx->methods,
                     encode->interpolation_methods[i],
                     0,
                     MAX_INTERPOLATIONS_EVER,
                     NULL
                );
        }
}
void input_interpolation_methods(coder_t *decode) {
        decode->nb_interpolation_methods
        = symb_get_int(&decode->coder,
                       &decode->ctx->methods,
                       0,
                       MAX_INTERPOLATIONS_EVER
                      ) + 1;
        if (decode->nb_interpolation_methods > MAX_INTERPOLATIONS) {
                fprintf(stderr,"WARNING: using more interpolation methods (%i) than supported by this implementation (%i).\n",decode->nb_interpolation_methods, MAX_INTERPOLATIONS);
        }
        for(int i=0; i<decode->nb_interpolation_methods; i++) {
                decode->interpolation_methods[i] 
                = symb_get_int(&decode->coder,
                     &decode->ctx->methods,
                     0,
                     MAX_INTERPOLATIONS_EVER
                );
                if (decode->interpolation_methods[i] > MAX_INTERPOLATIONS-1) {
                        fprintf(stderr,"WARNING: using interpolation method %i, which is not supported by this implementation.\n",decode->interpolation_methods[i]);
                }
        }
}
void output_interpolation_choice(coder_t *encode, int i, int c, int method, int prev) {
     int j=0;
     while( j<method ) {
        symb_put_bit(&encode->coder,&encode->ctx->interpolationMethod[i][c][j][prev],1,&encode->ctx->bitsInterpol[c]);
        encode->ctx->symbInterpol[c]++;
        j++;
     }
     if(j < encode->nb_interpolation_methods-1) {
       symb_put_bit(&encode->coder,&encode->ctx->interpolationMethod[i][c][j][prev],0,&encode->ctx->bitsInterpol[c]);
       encode->ctx->symbInterpol[c]++;
     }
}

int input_interpolation_choice(coder_t *decode, int i, int c, int prev) {
     int j=0;
     while ( j < decode->nb_interpolation_methods-1 
             && symb_get_bit(&decode->coder,&decode->ctx->interpolationMethod[i][c][j][prev])) {
        j++;
     }
     return j;
}

#define SPLIT_AVG_WINDOW 1
void update_avg_split_depths(coder_t *encode, int c, int depth, int line) {
        int *current_avg = (line? &encode->avg_Lsplit_depth[c] : &encode->avg_Asplit_depth[c]);
//        if (line) {fprintf(stderr,"Line "); } else {fprintf(stderr,"Area "); }
        
//        fprintf(stderr,"current avg split depth [%c] : %i \t(%i)\t",planesYIQ[c], *current_avg, depth);
        *current_avg = depth; //10*depth + rdiv(*current_avg*(SPLIT_AVG_WINDOW-10),SPLIT_AVG_WINDOW);
//        fprintf(stderr,"-> %i\n", *current_avg);
}

void output_mask(coder_t *encode,int oldmask,int newmask,int qsize,int depth,uint32_t *bvar, uint32_t *bdiffE, int line, const pixel_t *left, const pixel_t *right) {
     assert(oldmask <= newmask);
     int lrdiff = 0;
#if SPLIT_CHANNELS_SEPARATELY     
     for (int c = 0; c < 3 ; c++) {
#else
     int c = 0;
     {
#endif
       if (!(oldmask & (1 << c))) {
       int others=0;
       if (c==0) others = oldmask/2;
       if (c==1) others = ((oldmask & 4)>>1) | (newmask & 1);
       if (c==2) others = newmask & 3;
       if (others< 0 || others > 3) fprintf(stderr,"oei %i\n",others);
       if (line) { lrdiff = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]),(c==0 ? 255 : 510),0,QZ_LRDIFF_S); }
       else { lrdiff = quantize_log_uint32(bdiffE[c],QZ_LRDIFF_S); };
       int bit = 0;
       if ((oldmask & (1 << c)) != (newmask & (1 << c))) {
          bit = 1;
       } else {
          encode->splits[c][line]++;
       }
       encode->ctx->symbSplit[c][line]++;
//       fprintf(stderr,"BDIFF(%i): %i (stop:%i)\n",c,lrdiff,bit);

//        int current_avg = rdiv((line? encode->avg_Lsplit_depth[c] : encode->avg_Asplit_depth[c]),SPLIT_AVG_WINDOW);
#if SPLIT_CHANNELS_SEPARATELY     
        int properties[NB_PROPERTIES_S] = {qsize,lrdiff,qbvar(bvar[c],QZ_BVAR_S),others};
#else
//        int properties[NB_PROPERTIES_S] = {qsize,lrdiff,qbvar(bvar[c],QZ_BVAR_S),others};
        int properties[NB_PROPERTIES_S] = {qsize,qbvar(bvar[0],QZ_BVAR_S),qbvar((bvar[1]+bvar[2])/2,QZ_BVAR_S),others};
#endif
        int vcontext = NB_PROPERTIES_S;

        
#if FULL_BORDERS        
        context_leaf_bit *leaf = 
                        (line ? find_context_bit(encode->ctx->splitContL[c],properties,encode,1,c,&vcontext) 
                              : find_context_bit(encode->ctx->splitContA[c],properties,encode,2,c,&vcontext));
#else
        context_leaf_bit *leaf = find_context_bit(encode->ctx->splitContA[c],properties,encode,2,c,&vcontext);
#endif
//        symb_put_simple_bit(&encode->coder,bit,&encode->ctx->bitsSplit[c][line]);

        if (vcontext == NB_PROPERTIES_S) {
            uint64_t oldbits = encode->ctx->bitsSplit[c][line];
                symb_put_bit(&encode->coder,
                     &leaf->context,
                     bit,&encode->ctx->bitsSplit[c][line]);
            leaf->bits += encode->ctx->bitsSplit[c][line]  - oldbits;
        } else {
            uint64_t oldbits = leaf->splitbits[vcontext];
            int index=vcontext;
            if (properties[vcontext] > average(leaf->sum[vcontext],leaf->count)) index += NB_PROPERTIES_S;
            symb_put_bit(&encode->coder,
                     &leaf->splitcontext[index],
                     bit,
                     &leaf->splitbits[vcontext]);
            encode->ctx->bitsSplit[c][line] += leaf->splitbits[vcontext] - oldbits;
            symb_put_bit_v(&encode->coder,
                     &leaf->context,
                     bit,
                     &leaf->bits,
                     1 );
        }
        
        for (int i=0; i<NB_PROPERTIES_S ; i++) {
          if(leaf->splitbits[i]>0) {
           int index=i;
           if (properties[i] > average(leaf->sum[i],leaf->count)) index += NB_PROPERTIES_S;
           if (vcontext != i) {
             symb_put_bit_v(&encode->coder,
                     &leaf->splitcontext[index],
                     bit,
                     &leaf->splitbits[i],
                     1
                    );
           }
          }
        }

       
     }
    }
}

int input_mask(coder_t *decode,int mask,int qsize,int depth,uint32_t *bvar, uint32_t *bdiffE, int line, const pixel_t *left, const pixel_t *right) {
     int oldmask = mask;
     int newmask = oldmask;
     int lrdiff = 0;
#if SPLIT_CHANNELS_SEPARATELY     
     for (int c = 0; c < 3 ; c++) {
#else
     int c = 0;
     {
#endif
      if (!(oldmask & (1 << c))) {
       int others=0;
       if (c==0) others = oldmask/2;
       if (c==1) others = ((oldmask & 4)>>1) | (newmask & 1);
       if (c==2) others = newmask & 3;
       if (line) {lrdiff = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]),(c==0 ? 255 : 510),0,QZ_LRDIFF_S); }
       else { lrdiff = quantize_log_uint32(bdiffE[c],QZ_LRDIFF_S); };

#if SPLIT_CHANNELS_SEPARATELY     
        int properties[NB_PROPERTIES_S] = {qsize,lrdiff,qbvar(bvar[c],QZ_BVAR_S),others};
#else
//        int properties[NB_PROPERTIES_S] = {qsize,lrdiff,qbvar(bvar[c],QZ_BVAR_S),others};
        int properties[NB_PROPERTIES_S] = {qsize,qbvar(bvar[0],QZ_BVAR_S),qbvar((bvar[1]+bvar[2])/2,QZ_BVAR_S),others};
#endif
        int vcontext = NB_PROPERTIES_S;
        
#if FULL_BORDERS        
        context_leaf_bit *leaf= (line? find_context_bit(decode->ctx->splitContL[c],properties,decode,1,c,&vcontext)
                                 : find_context_bit(decode->ctx->splitContA[c],properties,decode,2,c,&vcontext));
#else
        context_leaf_bit *leaf= find_context_bit(decode->ctx->splitContA[c],properties,decode,2,c,&vcontext);
#endif
        int bit = 0;
        if (vcontext == NB_PROPERTIES_S) {
            uint64_t oldbits = decode->ctx->bitsSplit[c][line];
                bit = symb_get_bit_c(&decode->coder,
                     &leaf->context,
                     &decode->ctx->bitsSplit[c][line]);
            leaf->bits += decode->ctx->bitsSplit[c][line]  - oldbits;
        } else {
            uint64_t oldbits = leaf->splitbits[vcontext];
            int index=vcontext;
            if (properties[vcontext] > average(leaf->sum[vcontext],leaf->count)) index += NB_PROPERTIES_S;
            bit = symb_get_bit_c(&decode->coder,
                     &leaf->splitcontext[index],
                     &leaf->splitbits[vcontext]);
            decode->ctx->bitsSplit[c][line] += leaf->splitbits[vcontext] - oldbits;
            symb_put_bit_v(&decode->coder,
                     &leaf->context,
                     bit,
                     &leaf->bits,
                     1 );
        }
        for (int i=0; i<NB_PROPERTIES_S ; i++) {
          if(leaf->splitbits[i]>0) {
           int index=i;
           if (properties[i] > average(leaf->sum[i],leaf->count)) index += NB_PROPERTIES_S;
           if (vcontext != i) {
             symb_put_bit_v(&decode->coder,
                     &leaf->splitcontext[index],
                     bit,
                     &leaf->splitbits[i],
                     1
                    );
           }
          }
        }
        if (bit) newmask |= 1<<c;

       
     }
    }

#if SPLIT_CHANNELS_SEPARATELY == 0
    if (newmask>0) newmask=7;
#endif

    return newmask;
}

void output_min_max(coder_t *encode,int min,int max,int oldmin,int oldmax, int c) {
        assert(min <= max);
        symb_put_int(&encode->coder,
                     &encode->ctx->r_min[c],
                     min,
                     oldmin,
                     oldmax,
                     NULL
                    );
        symb_put_int(&encode->coder,
                     &encode->ctx->r_max[c],
                     (oldmax-max),
                     0,
                     (oldmax-min),
                     NULL
                    );
//        fprintf(stderr,"min: %i  (range: %i to %i)\n",min,oldmin,oldmax);                    
  //      fprintf(stderr,"max: %i  (range: %i to %i)\n",oldmax-max,0,oldmax-min);                    
    //    fprintf(stderr,"OUTPUT %i,%i,%i,%i,%i \n",min,max,oldmin,oldmax,c);                    
}
void input_min_max(coder_t *decode,int *min,int *max,int oldmin,int oldmax, int c) {
        *min = symb_get_int(&decode->coder,
                     &decode->ctx->r_min[c],
                     oldmin,
                     oldmax
                    );
        *max = oldmax - symb_get_int(&decode->coder,
                     &decode->ctx->r_max[c],
                     0,
                     (oldmax - *min)
                    );
//        fprintf(stderr,"min: %i   max: %i\n",*min,*max);                    
  //      fprintf(stderr,"min: %i  (range: %i to %i)\n",*min,oldmin,oldmax);                    
    //    fprintf(stderr,"max: %i  (range: %i to %i)\n",oldmax-*max,0,oldmax-*min);                    
      //  fprintf(stderr,"INPUT %i,%i,%i,%i,%i \n",*min,*max,oldmin,oldmax,c);                    
}



bitvector bv = {.nb = -1};


// output the color values (from planes not given by mask) of a pixel to the range encoder
void output_pixel(coder_t *encode, int x, int y, const pixel_t *guess, uint32_t size, int mask, 
                                const pixel_t *topleft, const pixel_t *topright, const pixel_t *left,const pixel_t *right, uint32_t bvar[3], uint32_t bdiffE[3], int depth, int onepixel, int splitdir) {
//  fprintf(stderr,"pixel (%i,%i)...\n",x,y);

//  pixels_at_depth[depth]++;
  const pixel_t def = {.d= {128,256,256}};
  if (!guess) guess = &def;

  int guesscolor[3];
  for (int c = 0; c<3; c++) guesscolor[c] = STRETCH(UNSTRETCH(guess->d[c]));
//  for (int c = 0; c<3; c++) guesscolor[c] = guess->d[c];


  pixel_t *p = image_pixel(encode->source,x,y);
  pixel_t *r = image_pixel(encode->rec,x,y);
  if (!topleft) topleft = &def;
  if (!topright) topright = &def;
  if (!left) left = &def;
  if (!right) right = &def;
  int qsize  = quantize_log_uint32(size,QZ_SIZE);
//  fprintf(stderr,"size = %i\n",size);
//  fprintf(stderr,"qsize = %i\n",qsize);
//  int qsize  = clamp(size-2,0,QZ_SIZE-1);
  for (int c = 0; c < 3; c++) {
      if (!(mask & (1 << c))) {
        encode->outputted_pixels[c]++;
#if (ACCURATE_YIQ_RANGES == 1)
        int minval=0;
        int maxval=0;

        switch (c) {
          case 0:
            get_range_y(&minval,&maxval,encode->range_min,encode->range_max);
            break;
          case 1:
            get_range_i(r->d[0],&minval,&maxval,encode->range_min,encode->range_max);
            break;
          case 2:
            get_range_q(r->d[0],r->d[1],&minval,&maxval,encode->range_min,encode->range_max);
            break;
        }
        if (p->d[c] < minval) {//fprintf(stderr,"Strange: %i < %i (channel %i)\n",p->d[c],minval,c);
                               p->d[c] = minval; }
        if (p->d[c] > maxval) {//fprintf(stderr,"Strange: %i > %i (channel %i)\n",p->d[c],maxval,c);
                                p->d[c] = maxval;}
#else
        int minval = STRETCH(encode->range_min[c]);
        int maxval = STRETCH(encode->range_max[c]);
#endif

        int variance  = quantize_log_uint32(
                        (c==0 ? 4000 : 1000) *
                        (4*UNSTRETCH(topleft->d[c])*UNSTRETCH(topleft->d[c])
                        +4*UNSTRETCH(topright->d[c])*UNSTRETCH(topright->d[c])
                        +4*UNSTRETCH(left->d[c])*UNSTRETCH(left->d[c])
                        +4*UNSTRETCH(right->d[c])*UNSTRETCH(right->d[c])
                -
                        (UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])+UNSTRETCH(left->d[c])+UNSTRETCH(right->d[c]))
                      * (UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])+UNSTRETCH(left->d[c])+UNSTRETCH(right->d[c]))
                        ),
                        QZ_VAR);
        assert(variance>=0);


        int qz = encode->qz[c];



        if (qz > COLOR_STRETCH) qz += COLOR_STRETCH*(depth-10)/EXTRA_QUANTIZATION_DEPTH 
                                    + COLOR_STRETCH*(variance/10 - 10)/EXTRA_QUANTIZATION_VAR;
        if (qz > MAX_QUANTIZATION*COLOR_STRETCH) qz = MAX_QUANTIZATION*COLOR_STRETCH;
        if (qz < COLOR_STRETCH) qz = COLOR_STRETCH;



        int dist = -rdiv(guesscolor[c] - (p->d[c]-(qz/2)),qz);
        int min = -rdiv(guesscolor[c] - (minval-(qz/2)),qz);
        int max = -rdiv(guesscolor[c] - (maxval-(qz/2)),qz);

//        fprintf(stderr,"Computing bitvector for dist=%i..%i..%i (channel %i)\n",min,dist,max,c);

        if (min==max || (encode->color_data.used && compute_range(&encode->color_data, min, max, c, r, qz, guesscolor[c],&bv) == 1)) {
//                if (dist != get_first_in_range(&bv,min,max,min)) fprintf(stderr,"Strange: expected that %i == %i\n",dist,get_first_in_range(&bv,min,max,min));
                assert(dist == get_first_in_range(&bv,min,max));
        } else {

        assert(check_range(&bv,dist,dist));

        int yg = UNSTRETCH(guesscolor[0])*QZ_YG/256;
        int ig = UNSTRETCH(guesscolor[1])*QZ_IG/511;
        int qg = UNSTRETCH(guesscolor[2])*QZ_QG/511;
        int ydiff = 0, idiff=0;
        switch (c) {
          case 0:
                break;
          case 1:
                ydiff = quantize_log(UNSTRETCH(guesscolor[0] - r->d[0]),255,1,QZ_YD);
                break;
          case 2:
                ydiff = quantize_log(UNSTRETCH(guesscolor[0] - r->d[0]),255,1,QZ_YD);
                idiff = quantize_log(UNSTRETCH(guesscolor[1] - r->d[1]),510,1,QZ_ID);
                break;
        }

        int nb_properties = (onepixel ? NB_PROPERTIES_1 : NB_PROPERTIES);
        int bits_used = 0;
        int vcontext = nb_properties;
        int properties[nb_properties];
        properties[0] = idiff;
        properties[1] = ydiff;
        properties[2] = yg;
        properties[3] = ig;
        properties[4] = qg;
        context_leaf *leaf;
#if FULL_BORDERS        
        if (onepixel) {
                properties[5]=quantize_log(UNSTRETCH(image_pixel(encode->rec,x-1,y-1)->d[c] - image_pixel(encode->rec,x+1,y+1)->d[c]),(c==0 ? 255 : 510),1,QZ_LRDIFF);
                properties[6]=quantize_log(UNSTRETCH(image_pixel(encode->rec,x+1,y-1)->d[c] - image_pixel(encode->rec,x-1,y+1)->d[c]),(c==0 ? 255 : 510),1,QZ_LRDIFF);
                properties[7]=quantize_log(UNSTRETCH(image_pixel(encode->rec,x,y-1)->d[c] - image_pixel(encode->rec,x,y+1)->d[c]),(c==0 ? 255 : 510),1,QZ_LRDIFF);
                properties[8]=quantize_log(UNSTRETCH(image_pixel(encode->rec,x-1,y)->d[c] - image_pixel(encode->rec,x+1,y)->d[c]),(c==0 ? 255 : 510),1,QZ_LRDIFF); 
                properties[9] = variance;
                properties[10] = qbvar(bvar[c],QZ_BVAR);
                leaf= find_context(encode->ctx->diff1[c],properties,encode,3,c,&vcontext);
        } else
#endif
        {
                int lrdiff = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        //        int grad   = quantize_log(UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])-UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]), (c==0 ? 510 : 1020),0,QZ_GRAD);
        //        int cross  = quantize_log(UNSTRETCH(topleft->d[c])+UNSTRETCH(right->d[c])-UNSTRETCH(left->d[c])-UNSTRETCH(topright->d[c]), (c==0 ? 510 : 1020),0,QZ_CROSS);
        //        assert(grad>=0);
        //        assert(cross>=0);
    
                assert(lrdiff>=0);
                assert(lrdiff<QZ_LRDIFF);
        //        assert(grad<QZ_GRAD);
        //        assert(cross<QZ_CROSS);
                assert(qsize>=0);
                properties[5] = qsize;
                properties[6] = lrdiff;
                properties[7] = variance;
                properties[8] = qbvar(bvar[c],QZ_BVAR);
                properties[9] = quantize_log_uint32(bdiffE[c],QZ_BDIFFE);
                properties[10] = splitdir;
//                fprintf(stderr,"qBVar=%i (%u) // qBDiffE=%i (%u)\n",properties[8],bvar[c],properties[9],bdiffE[c]);
        }
        if (encode->phase == 2) {
                symb_put_int_limited_a(&encode->coder,
                     find_context_p2(encode->ctx->diff_p2[c],properties,encode,0,c),
                     &dist,min,max,
                     &encode->ctx->bitsPixeldata[c][onepixel],
                     encode->sb[c],&bv
                    );
          
        } else {

          leaf= find_context(encode->ctx->diff[c],properties,encode,0,c,&vcontext);

          if (vcontext == nb_properties) {
            uint64_t oldbits = encode->ctx->bitsPixeldata[c][onepixel];
                symb_put_int_limited_a(&encode->coder,
                     &leaf->context,
                     &dist,min,max,
                     &encode->ctx->bitsPixeldata[c][onepixel],
                     encode->sb[c],&bv
//                     &encode->color_data, c, r, qz, guesscolor[c]
                    );
            bits_used = encode->ctx->bitsPixeldata[c][onepixel]  - oldbits;
            if(leaf->count >= CONTEXT_VIRTUAL_OUTPUT_DELAY) leaf->bits +=  bits_used;
          } else {
            uint64_t oldbits = leaf->splitbits[vcontext];
            int index=vcontext;
            if (properties[vcontext] > average(leaf->sum[vcontext],leaf->count)) index += nb_properties;
            symb_put_int_limited_a(&encode->coder,
                     &leaf->splitcontext[index],
                     &dist,min,max,
                     &leaf->splitbits[vcontext],
                     encode->sb[c],&bv
//                     &encode->color_data, c, r, qz, guesscolor[c]
                    );
            bits_used = leaf->splitbits[vcontext] - oldbits;
            encode->ctx->bitsPixeldata[c][onepixel] += bits_used;
            symb_put_int_limited_a1(&encode->coder,
                     &leaf->context,
                     &dist,min,max,
                     &leaf->bits,
                     encode->sb[c],
                     &bv
//                     &encode->color_data, c, r, qz, guesscolor[c]
                    );
          }
          if(leaf->count >= CONTEXT_VIRTUAL_OUTPUT_DELAY ) {
           for (int i=2-c; i<nb_properties ; i++) {
            if(leaf->splitbits[i]>0) {
             int index=i;
             if (properties[i] > average(leaf->sum[i],leaf->count)) index += nb_properties;
             if (vcontext != i) {
               symb_put_int_limited_a1(&encode->coder,
                     &leaf->splitcontext[index],
                     &dist,min,max,
                     &leaf->splitbits[i],
                     encode->sb[c],
                     &bv
//                     &encode->color_data, c, r, qz, guesscolor[c]
                    );
             }
            }
           }
          }
        }
        }
        
        
        encode->ctx->symbPixeldata[c][onepixel]++;
        r->d[c] = STRETCH(clamp(UNSTRETCH(dist*qz + guesscolor[c]),0,(c == 0 ? 255 : 510)));
//        if (p->d[c] != r->d[c]) fprintf(stderr,"CHANGE: %i became %i (channel %i)\n",  p->d[c],r->d[c],c);
        p->d[c] = r->d[c];
//        STRETCH(UNSTRETCH(dist*qz + guesscolor[c]));

//        p->d[c] = STRETCH((c==0?255:255));


// output bits per pixel
/*
        if (c==0) {
                p->d[0] = STRETCH(clamp(bits_used*bits_used/5461/5461,0,255));
        } else {
                if (mask&1) p->d[0]=0;
                p->d[0] = STRETCH(clamp(UNSTRETCH(p->d[0])+bits_used*bits_used/5461/5461,0,255));
                p->d[c] = STRETCH(256);
        }
*/


//        r->d[c] = dist*qz + guesscolor[c];


  //uncomment to see trajectory
/*
        if(depth==17){
//        if(qsize==0){
        r->d[0] = path[0];
        r->d[1] = STRETCH(150);
        r->d[2] = STRETCH(350);
        pixelcounter++;
        //if (pixelcounter%30 == 0) 
         {path[0] = (path[0]+1*COLOR_STRETCH)%32000;}
        }
*/
        
//        assert(encode->maxD[c]==0 ? (COLOR_DIST(r->d[c],p->d[c])==0) : 1);
      }
    }
    if (mask==0 && encode->qz[0] > COLOR_STRETCH) {
            round_color(&encode->color_data,r);
            *p = *r;
    }
#ifdef DEBUGMODE
//  fprintf(stderr,"%.*s",depth*2+1,"                                                                               ");
//  fprintf(stderr,"pixel (%i,%i) outputted\n",x,y);
#endif
}


#if FFV1
void output_pixel_ffv1(coder_t *encode, int x, int y, const pixel_t *guess, int mask,
                                const pixel_t *topleft, const pixel_t *left, const pixel_t *top,const pixel_t *leftleft,const pixel_t *toptop,const pixel_t *topright) {
  const pixel_t def = {.d= {128,256,256}};
  if (!guess) guess = &def;

  int guesscolor[3];
  for (int c = 0; c<3; c++) guesscolor[c] = STRETCH(UNSTRETCH(guess->d[c]));

  pixel_t *p = image_pixel(encode->source,x,y);
  pixel_t *r = image_pixel(encode->rec,x,y);
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        encode->outputted_pixels[c]++;
#if (ACCURATE_YIQ_RANGES == 1)
        int minval=0;
        int maxval=0;
        switch (c) {
          case 0:
            get_range_y(&minval,&maxval,encode->range_min,encode->range_max);
            break;
          case 1:
            get_range_i(r->d[0],&minval,&maxval,encode->range_min,encode->range_max);
            break;
          case 2:
            get_range_q(r->d[0],r->d[1],&minval,&maxval,encode->range_min,encode->range_max);
            break;
        }
        if (p->d[c] < minval) p->d[c] = minval;
        if (p->d[c] > maxval) p->d[c] = maxval;
#else
        int minval = STRETCH(encode->range_min[c]);
        int maxval = STRETCH(encode->range_max[c]);
#endif
        int qz = encode->qz[c];
        int dist = -rdiv(guesscolor[c] - (p->d[c]-(qz/2)),qz);
        int min = -rdiv(guesscolor[c] - (minval-(qz/2)),qz);
        int max = -rdiv(guesscolor[c] - (maxval-(qz/2)),qz);

        //bitvector bv = {};
//        int nb_possibilities = compute_range(&encode->color_data, min, max, c, r, qz, guesscolor[c],&bv);
//        if (nb_possibilities == 1) {
        if (min==max || (encode->color_data.used && compute_range(&encode->color_data, min, max, c, r, qz, guesscolor[c],&bv) == 1)) {
//                if (dist != get_first_in_range(&bv,min,max,min)) fprintf(stderr,"Strange: expected that %i == %i\n",dist,get_first_in_range(&bv,min,max,min));
                assert(dist == get_first_in_range(&bv,min,max));
        } else {
//        int yg = rdiv(guesscolor[0],qz)*QZ_YG/256;
//        int ig = rdiv(guesscolor[1],qz)*QZ_IG/511;
//        int qg = rdiv(guesscolor[2],qz)*QZ_QG/511;
        int yg = UNSTRETCH(guesscolor[0])*QZ_YG/256;
        int ig = UNSTRETCH(guesscolor[1])*QZ_IG/511;
        int qg = UNSTRETCH(guesscolor[2])*QZ_QG/511;
        int ydiff = 0, idiff=0;
        switch (c) {
          case 0:
                break;
          case 1:
                ydiff = quantize_log(UNSTRETCH(guesscolor[0] - r->d[0]),255,1,QZ_YD);
                break;
          case 2:
                ydiff = quantize_log(UNSTRETCH(guesscolor[0] - r->d[0]),255,1,QZ_YD);
                idiff = quantize_log(UNSTRETCH(guesscolor[1] - r->d[1]),510,1,QZ_ID);
                break;
        }
        int nb_properties = NB_PROPERTIES_FFV1;
        int bits_used = 0;
        int vcontext = nb_properties;
        int properties[nb_properties];
        properties[0] = idiff;
        properties[1] = ydiff;
        properties[2] = yg;
        properties[3] = ig;
        properties[4] = qg;
        properties[5] = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(topleft->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[6] = quantize_log(UNSTRETCH(topleft->d[c])-UNSTRETCH(top->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[7] = quantize_log(UNSTRETCH(top->d[c])-UNSTRETCH(topright->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[8] = quantize_log(UNSTRETCH(toptop->d[c])-UNSTRETCH(top->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[9] = quantize_log(UNSTRETCH(leftleft->d[c])-UNSTRETCH(left->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);

        if (encode->phase == 2) {
                symb_put_int_limited_a(&encode->coder,
                     find_context_p2(encode->ctx->diff_p2[c],properties,encode,0,c),
                     &dist,min,max,
                     &encode->ctx->bitsPixeldata[c][0],
                     encode->sb[c],&bv
                    );
          
        } else {



        context_leaf *leaf;
        leaf= find_context(encode->ctx->diff[c],properties,encode,0,c,&vcontext);

        if (vcontext == nb_properties) {
            uint64_t oldbits = encode->ctx->bitsPixeldata[c][0];
                symb_put_int_limited_a(&encode->coder,
                     &leaf->context,
                     &dist,min,max,
                     &encode->ctx->bitsPixeldata[c][0],
                     encode->sb[c],&bv
//                     &encode->color_data, c, r, qz, guesscolor[c]
                    );
            bits_used = encode->ctx->bitsPixeldata[c][0]  - oldbits;
            if(leaf->count >= CONTEXT_VIRTUAL_OUTPUT_DELAY) leaf->bits +=  bits_used;
        } else {
            uint64_t oldbits = leaf->splitbits[vcontext];
            int index=vcontext;
            if (properties[vcontext] > average(leaf->sum[vcontext],leaf->count)) index += nb_properties;
            symb_put_int_limited_a(&encode->coder,
                     &leaf->splitcontext[index],
                     &dist,min,max,
                     &leaf->splitbits[vcontext],
                     encode->sb[c],&bv
//                     &encode->color_data, c, r, qz, guesscolor[c]
                    );
            bits_used = leaf->splitbits[vcontext] - oldbits;
            encode->ctx->bitsPixeldata[c][0] += bits_used;
            symb_put_int_limited_a1(&encode->coder,
                     &leaf->context,
                     &dist,min,max,
                     &leaf->bits,
                     encode->sb[c],
                     &bv
//                     &encode->color_data, c, r, qz, guesscolor[c]
                    );
        }
        if(leaf->count >= CONTEXT_VIRTUAL_OUTPUT_DELAY) {
         for (int i=2-c; i<nb_properties ; i++) {
          if(leaf->splitbits[i]>0) {
           int index=i;
           if (properties[i] > average(leaf->sum[i],leaf->count)) index += nb_properties;
           if (vcontext != i) {
             symb_put_int_limited_a1(&encode->coder,
                     &leaf->splitcontext[index],
                     &dist,min,max,
                     &leaf->splitbits[i],
                     encode->sb[c],
                     &bv
//                     &encode->color_data, c, r, qz, guesscolor[c]
                    );
           }
          }
         }
        }
        }
        }
        encode->ctx->symbPixeldata[c][0]++;

//        r->d[c] = STRETCH(UNSTRETCH(dist*qz + guesscolor[c]));
        r->d[c] = STRETCH(clamp(UNSTRETCH(dist*qz + guesscolor[c]),0,(c == 0 ? 255 : 510)));
        p->d[c] = r->d[c];
/*        if (c==0) {
                p->d[0] = STRETCH(clamp(bits_used*bits_used/5461/5461,0,255));
        } else {
                if (mask&1) p->d[0]=0;
                p->d[0] = STRETCH(clamp(UNSTRETCH(p->d[0])+bits_used*bits_used/5461/5461,0,255));
                p->d[c] = STRETCH(256);
        }                
*/
        
      }
    }
//    round_color(&encode->color_data,r);
//    round_color(&encode->color_data,p);
    *p = *r;
#ifdef DEBUGMODE
//  fprintf(stderr,"%.*s",depth*2+1,"                                                                               ");
//  fprintf(stderr,"pixel (%i,%i) outputted: r=(%i,%i,%i), p=(%i,%i,%i)\n",x,y,r->d[0],r->d[1],r->d[2],p->d[0],p->d[1],p->d[2]);
#endif
//  fprintf(stderr,"pixel (%i,%i): r=(%i,%i,%i), guess=(%i,%i,%i)\n",x,y,r->d[0],r->d[1],r->d[2],guesscolor[0],guesscolor[1],guesscolor[2]);

}

void input_pixel_ffv1(coder_t *decode, int x, int y, const pixel_t *guess, int mask,
                                const pixel_t *topleft, const pixel_t *left, const pixel_t *top,const pixel_t *leftleft,const pixel_t *toptop,const pixel_t *topright) {
  const pixel_t def = {.d= {128,256,256}};
  if (!guess) guess = &def;

  int guesscolor[3];
  for (int c = 0; c<3; c++) guesscolor[c] = STRETCH(UNSTRETCH(guess->d[c]));

  pixel_t *r = image_pixel(decode->rec,x,y);
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
//        decode->outputted_pixels[c]++;
#if (ACCURATE_YIQ_RANGES == 1)
        int minval=0;
        int maxval=0;
        switch (c) {
          case 0:
            get_range_y(&minval,&maxval,decode->range_min,decode->range_max);
            break;
          case 1:
            get_range_i(r->d[0],&minval,&maxval,decode->range_min,decode->range_max);
            break;
          case 2:
            get_range_q(r->d[0],r->d[1],&minval,&maxval,decode->range_min,decode->range_max);
            break;
        }
#else
        int minval = STRETCH(decode->range_min[c]);
        int maxval = STRETCH(decode->range_max[c]);
#endif
        int qz = decode->qz[c];
        int dist = 0;
        int min = -rdiv(guesscolor[c] - (minval-(qz/2)),qz);
        int max = -rdiv(guesscolor[c] - (maxval-(qz/2)),qz);

        //bitvector bv = {};
//        int nb_possibilities = compute_range(&decode->color_data, min, max, c, r, qz, guesscolor[c],&bv);
//        if (nb_possibilities == 1) {
        if (min==max || (decode->color_data.used && compute_range(&decode->color_data, min, max, c, r, qz, guesscolor[c],&bv) == 1)) {
                dist = get_first_in_range(&bv,min,max);
        } else {

        int yg = UNSTRETCH(guesscolor[0])*QZ_YG/256;
        int ig = UNSTRETCH(guesscolor[1])*QZ_IG/511;
        int qg = UNSTRETCH(guesscolor[2])*QZ_QG/511;
      int ydiff = 0, idiff=0;
        switch (c) {
          case 0:
                break;
          case 1:
                ydiff = quantize_log(UNSTRETCH(guesscolor[0] - r->d[0]),255,1,QZ_YD);
                break;
          case 2:
                ydiff = quantize_log(UNSTRETCH(guesscolor[0] - r->d[0]),255,1,QZ_YD);
                idiff = quantize_log(UNSTRETCH(guesscolor[1] - r->d[1]),510,1,QZ_ID);
                break;
        }
        int nb_properties = NB_PROPERTIES_FFV1;
        int properties[nb_properties];
        properties[0] = idiff;
        properties[1] = ydiff;
        properties[2] = yg;
        properties[3] = ig;
        properties[4] = qg;
        properties[5] = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(topleft->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[6] = quantize_log(UNSTRETCH(topleft->d[c])-UNSTRETCH(top->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[7] = quantize_log(UNSTRETCH(top->d[c])-UNSTRETCH(topright->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[8] = quantize_log(UNSTRETCH(toptop->d[c])-UNSTRETCH(top->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[9] = quantize_log(UNSTRETCH(leftleft->d[c])-UNSTRETCH(left->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);

      if (decode->phase == 2) {
               dist= symb_get_int_limited_c_a(&decode->coder,
                     find_context_p2(decode->ctx->diff_p2[c],properties,decode,0,c),
                     min,max,
                     decode->sb[c],
                     &decode->ctx->bitsPixeldata[c][0],&bv                    );
          
      } else {


        int vcontext = nb_properties;
        context_leaf *leaf;
        leaf= find_context(decode->ctx->diff[c],properties,decode,0,c,&vcontext);

        if (vcontext == nb_properties) {
          uint64_t oldbits = decode->ctx->bitsPixeldata[c][0];
          dist = symb_get_int_limited_c_a(                &decode->coder,                &leaf->context,
                min,max,                decode->sb[c],                &decode->ctx->bitsPixeldata[c][0],&bv     );
          if(leaf->count >= CONTEXT_VIRTUAL_OUTPUT_DELAY) leaf->bits += decode->ctx->bitsPixeldata[c][0]  - oldbits;
        } else {
          uint64_t oldbits = leaf->splitbits[vcontext];
          int index=vcontext;
          if (properties[vcontext] > average(leaf->sum[vcontext],leaf->count)) index += nb_properties;
          dist = symb_get_int_limited_c_a(                &decode->coder,                &leaf->splitcontext[index],
                       min,max,                decode->sb[c],                &leaf->splitbits[vcontext],&bv       );
          decode->ctx->bitsPixeldata[c][0] += leaf->splitbits[vcontext]  - oldbits;
          symb_put_int_limited_a1(&decode->coder,                &leaf->context,
                &dist,min,max,                &leaf->bits,                decode->sb[c],                &bv        );
        }
        if(leaf->count >= CONTEXT_VIRTUAL_OUTPUT_DELAY) {
         for (int i=2-c; i<nb_properties ; i++) {
          if(leaf->splitbits[i]>0) {
           int index=i;
           if (properties[i] > average(leaf->sum[i],leaf->count)) index += nb_properties;
           if (vcontext != i) {
             symb_put_int_limited_a1(&decode->coder,                     &leaf->splitcontext[index],
                     &dist,min,max,       &leaf->splitbits[i],         decode->sb[c],                    &bv         );
           }
          }
         }
        }
        }
      }
//        r->d[c] = STRETCH(UNSTRETCH(dist*qz + guesscolor[c]));
        r->d[c] = STRETCH(clamp(UNSTRETCH(dist*qz + guesscolor[c]),0,(c == 0 ? 255 : 510)));
        
      }
    }

//    round_color(&decode->color_data,r);

//  fprintf(stderr,"pixel (%i,%i): r=(%i,%i,%i), guess=(%i,%i,%i)\n",x,y,r->d[0],r->d[1],r->d[2],guesscolor[0],guesscolor[1],guesscolor[2]);

#ifdef DEBUGMODE
//  fprintf(stderr,"%.*s",depth*2+1,"                                                                               ");
//  fprintf(stderr,"pixel (%i,%i) inputted\n",x,y);
#endif
}


#endif


// input the color values (from planes not given by mask) of a pixel from the range encoder
void input_pixel(coder_t *decode, int x, int y, const pixel_t *guess, uint32_t size, int mask,
                                const pixel_t *topleft, const pixel_t *topright, const pixel_t *left,const pixel_t *right, uint32_t bvar[3], uint32_t bdiffE[3], int depth, int onepixel, int splitdir) {
  const pixel_t def = {.d= {128,256,256}};
  pixel_t *r = image_pixel(decode->rec,x,y);
  if (!guess) guess = &def;

  int guesscolor[3];
  for (int c = 0; c<3; c++) guesscolor[c] = STRETCH(UNSTRETCH(guess->d[c]));
//  for (int c = 0; c<3; c++) guesscolor[c] = guess->d[c];
  int qsize  = quantize_log_uint32(size,QZ_SIZE);
//  int qsize  = clamp(size-2,0,QZ_SIZE-1);

  if (!topleft) topleft = &def;
  if (!topright) topright = &def;
  if (!left) left = &def;
  if (!right) right = &def;
    for (int c = 0; c < 3; c++) {
      if (!(mask & (1 << c))) {
#if (ACCURATE_YIQ_RANGES == 1)
        int minval=0, maxval=0;
        switch (c) {
          case 0:
            get_range_y(&minval,&maxval,decode->range_min,decode->range_max);
            break;
          case 1:
            get_range_i(r->d[0],&minval,&maxval,decode->range_min,decode->range_max);
            break;
          case 2:
            get_range_q(r->d[0],r->d[1],&minval,&maxval,decode->range_min,decode->range_max);
            break;
        }
#else
        int minval = STRETCH(decode->range_min[c]);
        int maxval = STRETCH(decode->range_max[c]);
#endif
        int variance  = quantize_log_uint32(
                        (c==0 ? 4000 : 1000) *
                        (4*UNSTRETCH(topleft->d[c])*UNSTRETCH(topleft->d[c])
                        +4*UNSTRETCH(topright->d[c])*UNSTRETCH(topright->d[c])
                        +4*UNSTRETCH(left->d[c])*UNSTRETCH(left->d[c])
                        +4*UNSTRETCH(right->d[c])*UNSTRETCH(right->d[c])
                -
                        (UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])+UNSTRETCH(left->d[c])+UNSTRETCH(right->d[c]))
                      * (UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])+UNSTRETCH(left->d[c])+UNSTRETCH(right->d[c]))
                        ),
                        QZ_VAR);
        assert(variance>=0);


        int qz = decode->qz[c];
        if (qz > COLOR_STRETCH) qz += COLOR_STRETCH*(depth-10)/EXTRA_QUANTIZATION_DEPTH 
                                    + COLOR_STRETCH*(variance/10 - 10)/EXTRA_QUANTIZATION_VAR;
        if (qz > MAX_QUANTIZATION*COLOR_STRETCH) qz = MAX_QUANTIZATION*COLOR_STRETCH;
        if (qz < COLOR_STRETCH) qz = COLOR_STRETCH;



        int min = -rdiv(guesscolor[c] - (minval-(qz/2)),qz);
        int max = -rdiv(guesscolor[c] - (maxval-(qz/2)),qz);

        //bitvector bv = {};
        int dist=0;
        
//        int nb_possibilities = compute_range(&decode->color_data, min, max, c, r, qz, guesscolor[c],&bv);
//        if (nb_possibilities == 1) {
        if (min==max || (decode->color_data.used && compute_range(&decode->color_data, min, max, c, r, qz, guesscolor[c],&bv) == 1)) {
                dist = get_first_in_range(&bv,min,max);
        } else {
        assert(check_range(&bv,min,max));
        int yg = UNSTRETCH(guesscolor[0])*QZ_YG/256;
        int ig = UNSTRETCH(guesscolor[1])*QZ_IG/511;
        int qg = UNSTRETCH(guesscolor[2])*QZ_QG/511;
        int ydiff = 0, idiff=0;
        switch (c) {
          case 0:
                break;
          case 1:
                ydiff = quantize_log(UNSTRETCH(guesscolor[0] - r->d[0]),255,1,QZ_YD);
                break;
          case 2:
                ydiff = quantize_log(UNSTRETCH(guesscolor[0] - r->d[0]),255,1,QZ_YD);
                idiff = quantize_log(UNSTRETCH(guesscolor[1] - r->d[1]),510,1,QZ_ID);
                break;
        }

        int nb_properties = (onepixel ? NB_PROPERTIES_1 : NB_PROPERTIES);
        int vcontext = nb_properties;
        int properties[nb_properties];
        properties[0] = idiff;
        properties[1] = ydiff;
        properties[2] = yg;
        properties[3] = ig;
        properties[4] = qg;
        context_leaf *leaf;
#if FULL_BORDERS        
        if (onepixel) {
                properties[5]=quantize_log(UNSTRETCH(image_pixel(decode->rec,x-1,y-1)->d[c] - image_pixel(decode->rec,x+1,y+1)->d[c]),(c==0 ? 255 : 510),1,QZ_LRDIFF);
                properties[6]=quantize_log(UNSTRETCH(image_pixel(decode->rec,x+1,y-1)->d[c] - image_pixel(decode->rec,x-1,y+1)->d[c]),(c==0 ? 255 : 510),1,QZ_LRDIFF);
                properties[7]=quantize_log(UNSTRETCH(image_pixel(decode->rec,x,y-1)->d[c] - image_pixel(decode->rec,x,y+1)->d[c]),(c==0 ? 255 : 510),1,QZ_LRDIFF);
                properties[8]=quantize_log(UNSTRETCH(image_pixel(decode->rec,x-1,y)->d[c] - image_pixel(decode->rec,x+1,y)->d[c]),(c==0 ? 255 : 510),1,QZ_LRDIFF);
                properties[9] = variance;
                properties[10] = qbvar(bvar[c],QZ_BVAR);
                leaf= find_context(decode->ctx->diff1[c],properties,decode,3,c,&vcontext);
        } else 
#endif
        {
                int lrdiff = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        //        int grad   = quantize_log(UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])-UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]), (c==0 ? 510 : 1020),0,QZ_GRAD);
        //        int cross  = quantize_log(UNSTRETCH(topleft->d[c])+UNSTRETCH(right->d[c])-UNSTRETCH(left->d[c])-UNSTRETCH(topright->d[c]), (c==0 ? 510 : 1020),0,QZ_CROSS);
        //        assert(grad>=0);
        //        assert(cross>=0);

                assert(lrdiff>=0);
                assert(lrdiff<QZ_LRDIFF);
        //        assert(grad<QZ_GRAD);
        //        assert(cross<QZ_CROSS);
                assert(qsize>=0);
                properties[5] = qsize;
//              properties[6] = depth;
                properties[6] = lrdiff;
                properties[7] = variance;
                properties[8] = qbvar(bvar[c],QZ_BVAR);
                properties[9] = quantize_log_uint32(bdiffE[c],QZ_BDIFFE);
                properties[10] = splitdir;
//                properties[9] = x;
//                properties[10] = y;
        }
      if (decode->phase == 2) {
               dist= symb_get_int_limited_c_a(&decode->coder,
                     find_context_p2(decode->ctx->diff_p2[c],properties,decode,0,c),
                     min,max,
                     decode->sb[c],
                     &decode->ctx->bitsPixeldata[c][onepixel],&bv                    );
          
      } else {
        leaf= find_context(decode->ctx->diff[c],properties,decode,0,c,&vcontext);
        if (vcontext == nb_properties) {
          uint64_t oldbits = decode->ctx->bitsPixeldata[c][onepixel];
          dist = symb_get_int_limited_c_a(
                &decode->coder,                &leaf->context,                min,max,                decode->sb[c],
                &decode->ctx->bitsPixeldata[c][onepixel],&bv               );
          if(leaf->count >= CONTEXT_VIRTUAL_OUTPUT_DELAY) leaf->bits += decode->ctx->bitsPixeldata[c][onepixel]  - oldbits;
        } else {
          uint64_t oldbits = leaf->splitbits[vcontext];
          int index=vcontext;
          if (properties[vcontext] > average(leaf->sum[vcontext],leaf->count)) index += nb_properties;
          dist = symb_get_int_limited_c_a(
                &decode->coder,                &leaf->splitcontext[index],                min,max,                decode->sb[c],
                &leaf->splitbits[vcontext],&bv               );
          decode->ctx->bitsPixeldata[c][onepixel] += leaf->splitbits[vcontext]  - oldbits;
          symb_put_int_limited_a1(&decode->coder,
                &leaf->context,                &dist,min,max,                &leaf->bits,                decode->sb[c],
                &bv               );
        }
        if(leaf->count >= CONTEXT_VIRTUAL_OUTPUT_DELAY && decode->ctx->treesize[0][c] < MAX_CONTEXTS) {
         for (int i=2-c; i<nb_properties ; i++) {
          if(leaf->splitbits[i]>0) {
           int index=i;
           if (properties[i] > average(leaf->sum[i],leaf->count)) index += nb_properties;
           if (vcontext != i) {
             symb_put_int_limited_a1(&decode->coder,
                     &leaf->splitcontext[index],                     &dist,min,max,                     &leaf->splitbits[i],
                     decode->sb[c],                     &bv                    );
           }
          }
         }
        }
        }
        }
        r->d[c] = STRETCH(clamp(UNSTRETCH(dist*qz + guesscolor[c]),0,(c == 0 ? 255 : 510)));
//        r->d[c] = STRETCH(UNSTRETCH(dist*qz + guesscolor[c]));
      }
    }
    if (mask==0 && decode->qz[0] > COLOR_STRETCH) {
       round_color(&decode->color_data,r);
    }
}


void output_colors(coder_t *coder) {


    int fullbuckets = 0;
    int emptybuckets = 0;
    int partialbuckets = 0;
    int indexedcolors = 0;
    
    for (int i=coder->range_min[0] ; i <= coder->range_max[0] ; i++) {
//       symb_put_simple_int(&coder->coder,coder->color_data.bucketY[i],0,1,&(coder->ctx->bitsColorData));
       symb_put_int(&coder->coder,&coder->ctx->colors_Y,coder->color_data.bucketY[i],0,1,&(coder->ctx->bitsColorData_Y));
    }
    for (int i=coder->range_min[0] / SIZE_Y; i <= coder->range_max[0] / SIZE_Y; i++) {
      int min_y = i*SIZE_Y;
      for ( ; min_y < (i+1)*SIZE_Y ; min_y++ ) {
        if (coder->color_data.bucketY[min_y]) break;
      }
      int max_y = (i+1)*SIZE_Y-1;
      for ( ; max_y >= i*SIZE_Y ; max_y-- ) {
        if (coder->color_data.bucketY[max_y]) break;
      }
      if (min_y < (i+1)*SIZE_Y) {
       for (int j=coder->range_min[1] / SIZE_I; j <= coder->range_max[1] / SIZE_I; j++) {
        for (int k=coder->range_min[2] / SIZE_Q; k <= coder->range_max[2] / SIZE_Q; k++) {
          if (bucket_exists(i,j,k)) {
           int c = coder->color_data.bucketYIQ[i][j][k].count;
           if (c == 0) { emptybuckets++; }
           else if (c < MAX_PER_BUCKET) {partialbuckets++; indexedcolors += c;}

//           fprintf(stderr,"Bucket[%i,%i,%i].count = %i\n",i,j,k,c);
//           symb_put_simple_int_0(&coder->coder,c,0,MAX_PER_BUCKET,&(coder->ctx->bitsColorData));
           symb_put_int(&coder->coder,&coder->ctx->colors_count0,(c==0),0,1,&(coder->ctx->bitsColorData_count));
           if (c>0) {
                if (c==MAX_PER_BUCKET) {
                  symb_put_int(&coder->coder,&coder->ctx->colors_count,0,0,MAX_PER_BUCKET-1,&(coder->ctx->bitsColorData_count));
                } else {
                  symb_put_int(&coder->coder,&coder->ctx->colors_count,c,0,MAX_PER_BUCKET-1,&(coder->ctx->bitsColorData_count));
                }
           }
           int bb_ranges = 0;
           int prev_y = min_y;
           int prev_i = -1;
           int prev_q = -1;
           if (c == MAX_PER_BUCKET) {c=2; bb_ranges=1; fullbuckets++;}
           for (int n=0; n < c; n++) {
            int cy = coder->color_data.bucketYIQ[i][j][k].color[n].d[0];
            int ci = coder->color_data.bucketYIQ[i][j][k].color[n].d[1];
            int cq = coder->color_data.bucketYIQ[i][j][k].color[n].d[2];
            int lower_y = clamp(prev_y,i*SIZE_Y,(i+1)*SIZE_Y-1);
            symb_put_int(&coder->coder,&coder->ctx->colors[0],cy-lower_y,0,max_y-lower_y,&(coder->ctx->bitsColorData[0]));

            int lower_i = clamp((cy == prev_y ? prev_i : j*SIZE_I),j*SIZE_I,(j+1)*SIZE_I-1);
            if (bb_ranges && n==1) lower_i = prev_i;
            symb_put_int(&coder->coder,&coder->ctx->colors[1],ci-lower_i,0,(j+1)*SIZE_I-1-lower_i,&(coder->ctx->bitsColorData[1]));

            int lower_q = clamp(( (cy==prev_y && ci==prev_i) ? prev_q+1 : k*SIZE_Q),k*SIZE_Q,(k+1)*SIZE_Q-1);
            if (bb_ranges && n==1) lower_q = prev_q;
            symb_put_int(&coder->coder,&coder->ctx->colors[2],cq-lower_q,0,(k+1)*SIZE_Q-1-lower_q,&(coder->ctx->bitsColorData[2]));
            prev_y = cy;
            prev_i = ci;
            prev_q = cq;

//            fprintf(stderr,"Bucket[%i,%i,%i].color[%i] = (%i,%i,%i)\n",i,j,k,n,coder->color_data.bucketYIQ[i][j][k].color[n].d[0],coder->color_data.bucketYIQ[i][j][k].color[n].d[1],coder->color_data.bucketYIQ[i][j][k].color[n].d[2]);
//            fprintf(stderr,"Bucket[%i,%i].color[%i] = (%i,%i,%i)\n",i,j,n,coder->color_data.bucketYI[i][j].color[n].d[0],coder->color_data.bucketYI[i][j].color[n].d[1],coder->color_data.bucketYI[i][j].color[n].d[2]);
           }
          }
        }
       }
      }
    }
    fill_bucketYI(&coder->color_data,coder->range_min[2],coder->range_max[2]);
/*    fprintf(stderr,"Outputted %i buckets: %i empty, %i full, %i partial (containing %i exact indexed colors)\n",
     fullbuckets + emptybuckets +  partialbuckets,
     emptybuckets, fullbuckets, partialbuckets, indexedcolors);*/
     
}
void input_colors(coder_t *coder) {
    const pixel_t init_max = {.d={255,510,510}};
    for (int i=coder->range_min[0] ; i <= coder->range_max[0] ; i++) {
       coder->color_data.bucketY[i] = symb_get_int(&coder->coder,&coder->ctx->colors_Y,0,1);
    }
    for (int i=coder->range_min[0] / SIZE_Y; i <= coder->range_max[0] / SIZE_Y; i++) {
      int min_y = i*SIZE_Y;
      for ( ; min_y < (i+1)*SIZE_Y ; min_y++ ) {
        if (coder->color_data.bucketY[min_y]) break;
      }
      int max_y = (i+1)*SIZE_Y-1;
      for ( ; max_y >= i*SIZE_Y ; max_y-- ) {
        if (coder->color_data.bucketY[max_y]) break;
      }
      if (min_y < (i+1)*SIZE_Y) {
       for (int j=coder->range_min[1] / SIZE_I; j <= coder->range_max[1] / SIZE_I; j++) {
        for (int k=coder->range_min[2] / SIZE_Q; k <= coder->range_max[2] / SIZE_Q; k++) {
          if (bucket_exists(i,j,k)) {
//            int c = symb_get_simple_int_0(&coder->coder,0,MAX_PER_BUCKET);
            int c = symb_get_int(&coder->coder,&coder->ctx->colors_count0,0,1);
            if (c == 0) {
                c = symb_get_int(&coder->coder,&coder->ctx->colors_count,0,MAX_PER_BUCKET-1);
                if (c == 0) c=MAX_PER_BUCKET;
            } else {
                c = 0;
            }

//            fprintf(stderr,"Bucket[%i,%i,%i].count = %i\n",i,j,k,c);
            int bb_ranges = 0;
            if (c == 0 || c == MAX_PER_BUCKET) {
             coder->color_data.bucketYIQ[i][j][k].count = c;
             if (c == MAX_PER_BUCKET) {
               coder->color_data.bucketYIQ[i][j][k].color[0] = init_max;  // make sure lower bound gets adjusted when adding colors later
               c = 2;
               bb_ranges=1;
             }
            }

            int prev_y = min_y;
            int prev_i = -1;
            int prev_q = -1;
            for (int n=0; n < c; n++) {
             int lower_y = clamp(prev_y,i*SIZE_Y,(i+1)*SIZE_Y-1);
             int cy = symb_get_int(&coder->coder,&coder->ctx->colors[0],0,max_y-lower_y)+lower_y;

             int lower_i = clamp((cy == prev_y ? prev_i : j*SIZE_I),j*SIZE_I,(j+1)*SIZE_I-1);
             if (bb_ranges && n==1) lower_i = prev_i;
             int ci = symb_get_int(&coder->coder,&coder->ctx->colors[1],0,(j+1)*SIZE_I-1-lower_i)+lower_i;

             int lower_q = clamp(( (cy==prev_y && ci==prev_i) ? prev_q+1 : k*SIZE_Q),k*SIZE_Q,(k+1)*SIZE_Q-1);
             if (bb_ranges && n==1) lower_q = prev_q;
             int cq = symb_get_int(&coder->coder,&coder->ctx->colors[2],0,(k+1)*SIZE_Q-1-lower_q)+lower_q;

             pixel_t color = {.d={STRETCH(cy),STRETCH(ci),STRETCH(cq)}};
             add_color(&coder->color_data,&color);
             prev_y = cy;
             prev_i = ci;
             prev_q = cq;
//            fprintf(stderr,"Bucket[%i,%i,%i].color[%i] = (%i,%i,%i)\n",i,j,k,n,coder->color_data.bucketYIQ[i][j][k].color[n].d[0],coder->color_data.bucketYIQ[i][j][k].color[n].d[1],coder->color_data.bucketYIQ[i][j][k].color[n].d[2]);
//            fprintf(stderr,"Bucket[%i,%i].color[%i] = (%i,%i,%i)\n",i,j,n,coder->color_data.bucketYI[i][j].color[n].d[0],coder->color_data.bucketYI[i][j].color[n].d[1],coder->color_data.bucketYI[i][j].color[n].d[2]);
            }
          }
        }
       }
      }
    }
    fill_bucketYI(&coder->color_data,coder->range_min[2],coder->range_max[2]);
}

