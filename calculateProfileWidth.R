calcProfileWidth <- function(alpha, theta, r){
        if(alpha > 2*pi | alpha < 0) 
		stop('alpha is out of bounds. alpha should be in 0<a<2*pi')
        if(theta > 2*pi | theta < 0) 
		stop('theta is out of bounds. theta should be in 0<a<2*pi')

	if(alpha > pi){
	        if(alpha < 4*pi - 2*theta){
		        p <- r*(theta - cos(alpha/2) + 1)/pi
                } else if(alpha <= 3*pi - theta){
                        p <- r*(theta - cos(alpha/2) + cos(alpha/2 + theta))/pi
                } else {
                        p <- r*(theta + 2*sin(theta/2))/pi
                }
        } else {
        	if(alpha < 4*pi - 2*theta){
                        p <- r*(theta*sin(alpha/2) - cos(alpha/2) + 1)/pi
 		} else {
                        p <- r*(theta*sin(alpha/2) - cos(alpha/2) + cos(alpha/2 + theta))/pi
                }
        }
        return(p)
}