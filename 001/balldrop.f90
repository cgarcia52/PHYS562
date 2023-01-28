program ball

    !use numtype
    use drop_setup
    implicit none


    real(dp) :: t, dt, x, y, v0, height, alpha, alpha_rad, vx0, vy0, & 
    vx, vy

    namelist /ball_input/ height, alpha, v0
    open(unit = 10, file = 'ball_input.nml')
    read(10, nml = ball_input)
    close(unit = 10)

print *, height, v0, alpha

alpha_rad = alpha * pi/180

vx0 = v0 * cos(alpha_rad)
vy0 = v0 * sin(alpha_rad)

x = 0
y = height
t = 0 
dt = 0.01_dp

do while (y > 0)

    vx = vx0
    vy = vy0 - g*t

    x = vx0 * t 
    y = height + vy0* t - g/2 * t**2

    t = t + dt

    write(3, *) t, vx, vy
    write(4, *) t, x, y
    write(7, *) x, y 

end do


end program ball