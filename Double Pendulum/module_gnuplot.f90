MODULE gnuplot
             CONTAINS
  
             SUBROUTINE draw_pendulum(x1,x2,y1,y2,l1,l2,t,k)

             IMPLICIT NONE
             REAL*8 , INTENT(IN)   :: x1, x2, y1, y2, l1, l2, t        
             INTEGER, INTENT(IN)   :: k                                             
             REAL*8                :: l                                            
             CHARACTER(LEN=16)     :: fname1                                      
             CHARACTER(LEN=6)      :: string1                                      
             CHARACTER(LEN=16)     :: fname2                                       
             CHARACTER(LEN=10)     :: string2                                     



      l = l1 + l2 + 1.d0

      WRITE(string1,'(i6.6)') k                                                    
              fname1 = 'pend'//string1//'.png'                                       
      
      WRITE(string2,'(F10.6)') t                                        
              fname2 = string2                                                     
        

      OPEN (20,FILE = "script_pendulum.gnu")                                      
      
      WRITE(20,FMT = '(a)') 'reset'                                                   
      
      WRITE(20,FMT = '(a)') "set terminal pngcairo dashed size 700,700 enhanced font 'Verdana,20'"                                      
 
      WRITE(20,FMT='(a , a , a)')    "set output '",fname1,"'"                        

      WRITE(20,FMT = '(a)') ''                                                        

      WRITE(20,FMT = '(a)') 'unset key'                                                

      WRITE(20,FMT='(a , a , a)')    "set title ' Tempo = ",fname2," s'"              

      WRITE(20,FMT = '(a)') '' 

      WRITE(20,FMT = '(a)') ' set format x'' '                                       
 
      WRITE(20,FMT = '(a)') ' set format y'' '                                        

      WRITE(20,FMT = '(a)') '' 
 
      WRITE(20,FMT = '(a)') 'set size ratio -1'                                       

      WRITE(20,FMT = '(a)') 'set style fill solid 1.0 border -1'                      

      WRITE(20,FMT = '(a)') '' 


      WRITE(20,FMT = '(a,F4.2,a,F4.2,a)') 'set xrange[-',l,':',l,']'                  
      WRITE(20,FMT = '(a,F4.2,a,F4.2,a)') 'set yrange[-',l,':',l,']'                  


      WRITE(20,FMT = '(a)') '' 


      WRITE(20,FMT = '(a,F12.5,a,F12.5,a)') 'set arrow 1 from 0.0 , 0.0 to ', x1 ,',', y1 ,' nohead front lw 2'                                       
      
      WRITE(20,FMT = '(a,F12.5,a,F12.5,a,F12.5,a,F12.5,a)') 'set arrow 2 from ', x1 ,',', y1 ,' to ', x2 ,',', y2 ,' nohead front lw 2'                
      
      WRITE(20,FMT = '(a,F4.2,a,F4.2,a)') 'set arrow 3 from -',l,' , 0.0 to ', l ,',0.0 nohead front lw 2'                                             

      WRITE(20,FMT = '(a)') '' 




      WRITE(20,FMT = '(a,F12.5,a,F12.5,a)') 'set object 1  circle center ', x1 ,',', y1 ,' size 0.1 fc rgb "black" front lw 3'                        
      WRITE(20,FMT = '(a,F12.5,a,F12.5,a)') 'set object 2  circle center ', x2 ,',', y2 ,' size 0.1 fc rgb "black" front lw 3'                     



      WRITE(20,FMT = '(a)') ' '


      WRITE(20,FMT = '(a)') 'plot NaN'                                              
      WRITE(20,FMT = '(a)') 'reset'

      WRITE(20,FMT = '(a)') '' 
      WRITE(20,FMT = '(a)') '' 




             END SUBROUTINE draw_pendulum 
             
END MODULE gnuplot