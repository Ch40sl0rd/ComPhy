#include<stdio.h>
#include<stdlib.h>



typedef struct listnode
{
  struct listnode *next;
  struct listnode *prev;
  double data;
} listnode;

typedef struct list
{
  struct listnode *first;
  struct listnode *last;
} list;

list *listcreate(void){
   list *new=(list*)malloc(sizeof(list)*1);
   if (new==NULL){
     printf("Leider konnen wir nicht mehr speicher allokieren\n");
     exit(1);
   }
   new->first=NULL;
   new->last=NULL;
   return new;
}

listnode *list_insert(list *l, listnode *e, double p){
  int s=0;
  if (l==NULL){
    printf("Erstmal musstest du ein list erzeugen\n");
    exit(1);
  }
  listnode *new=(listnode *)malloc(sizeof(listnode)*1);
  if (new == NULL){
    printf("Es gibt keinen speicher mehr zur verfugung\n");
    exit(1);
  }
  new->data=p;
  new->prev=NULL;
  new->next=NULL;
  if (e==NULL){
    if (l->first == NULL){
      l->first= new;
      l->last = new;
    }
    else{
      new->next=l->first;
      l->first->prev=new;
      l->first=new;
    }
    return new;
  }
  if (e->next == NULL){
    e->next=new;
    l->last=new;
    return new;
  }
  if (e==l->last)
    s=1;
  listnode *temp=e->next;
  e->next=new;
  new->prev=e;
  new->next=temp;
  if (s==1){
   l->last=new;
  }  
  return new; 
}
listnode *list_unshift(list *l,double p){
  listnode *new=(listnode *)malloc(sizeof(listnode)*1);
  if (new == NULL){
    printf("Erstmal musstest du ein list erzeugen\n");
    exit(1);
  }
  new->data=p;
  new->next=l->first;
  new->prev=NULL;
  l->first=new;
  return new;
}

listnode *list_push(list *l,double p){
  listnode *new=(listnode *)malloc(sizeof(listnode)*1);
  if (new == NULL){
    printf("Erstmal musstest du ein list erzeugen\n");
    exit(1);
  }
  new->data=p;
  new->prev=l->last;
  new->next=NULL;
  l->last=new;
  return new;
}

double list_shift(list *l){
  double ret;
  if (l->first ==NULL){
    printf("Den list ist leer\n");
    exit(1);
  }
  listnode *temp=l->first;
  if (temp->next == NULL){
    ret=l->first->data;
    free(temp);
    l->first=NULL;
    l->last=NULL;
    return ret;
  }
  else{
    temp->next->prev=NULL;
    ret=l->first->data;
    l->first=temp->next;
    free(temp);
    return ret;
  }
}
double list_pop(list *l){
  double ret;
  if (l->last==NULL){
    printf("Den list ist leer\n");
    return -1;
  }
  listnode *temp=l->last;
  if (temp->prev == NULL){
    ret=l->last->data;
    free(temp);
    l->first=NULL;
    l->last=NULL;
    return ret;
  }
  else{
    temp->prev->next=NULL;
    ret=l->last->data;
    l->last=temp->prev;
    free(temp);
    return ret;
  }
}
void listprintbackward(list *l){
  listnode *temp=l->last;
  if (l->last ==NULL){
   printf("Den list ist leer\n");
   return;
  }
  do{
   printf("Element von list ist %e\n",temp->data);
   temp=temp->prev;
  }
  while( temp !=NULL);;
}

void listprint(list *l){
  listnode *temp=l->first;
  if (l->first ==NULL){
   printf("Den list ist leer\n");
   return;
  }
  do{
   printf("Element von list ist %e\n",temp->data);
   temp=temp->next;
  }
  while( temp !=NULL);;
}
void list_delete(list *l, listnode *e){

  if (e==l->first){
     if (e==l->last){
       free(e);
       l->first=NULL;
       l->last=NULL;
     }
     l->first=e->next;
     free(e);
     l->first->prev=NULL;
     return;
  }
  if (e==l->last){
     l->last=e->prev;
     free(e);
     l->last->next=NULL;
     return;
  }
  listnode *temp1=e->next;
  listnode *temp2=e->prev;
  e->prev->next=temp1;
  e->next->prev=temp2;
  free(e);
}
list *merge(list *a,list *b){
  list *new=listcreate();
  new->first=a->first;
  a->last->next=b->first;
  b->first->prev=a->last;
  new->last=b->last;
  return new;
}
void list_free(list *l){
  if (l->first == NULL){
    free(l);
    return;
  }
  listnode *temp=l->first;
  listnode *temp2;
  do{
   temp2=temp->next;
   free(temp);
   temp=temp2;
  }
  while( temp!=NULL );
  free(l);
  return;
}
int main(int argc, char *argv[]){

   list *test;
   test=listcreate();
   list_insert(test, NULL, 10.);
//   if (test->first->prev==NULL){
//       printf("Ennek igy jell lennie\n");
//      }

   list_insert(test, NULL, 12.);
   list_insert(test, NULL, 13.);
   listprint(test);

   printf("Removing from the end\n");
   double tmp=list_pop(test);
   printf("The element that were removed %e\n",tmp);
   printf("The remaining list\n");
   //listprintbackward(test);
   printf("Leerzeile\n");
   list_free(test);
   listprint(test);
}
