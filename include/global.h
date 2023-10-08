//Includes Global Variables
#pragma once
class Parameters
{
    public:
        static Parameters& getInstance()
        {
            static Parameters    instance; // Guaranteed to be destroyed.
                                  // Instantiated on first use.
            return instance;
        }
    private:
        Parameters() {}            

    public:
        Parameters(Parameters const&)      = delete;
        void operator=(Parameters const&)  = delete;
        
        bool logData = false;
        unsigned int mode = 0;
        bool stereo = 0;
        float rf_fs, if_fs, audio_fs;
        bool rds = false;
        int sps = 0;
        int rf_down = 0, rf_up =0,if_up = 0,if_down=0;

        int stereo_lower_freq, stereo_upper_freq=0;
	    float carrier_lower_freq=0, carrier_upper_freq=0;

        int rdsUp = 0;
        int rdsDown = 0;

        int blockSize = 0;
};