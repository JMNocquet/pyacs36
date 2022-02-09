###################################################################
def save_velocity(self,vel_file='../stat/vel'):
###################################################################
    """
    save horizontal and up velocities
    """
    h_vel=vel_file+'_en.gmt'
    u_vel=vel_file+'_up.gmt'
    for gts in self.lGts():
        try:
            gts.save_velocity(h_vel)
        except:
            print("!!! Could not save ",gts.code," in ",h_vel)
        try:
            gts.save_velocity(u_vel,up=True)
        except:
            print("!!! Could not save ",gts.code," in ",u_vel)
    return()
