# Coral Bleaching Dataset --- Data Dictionary

This data dictionary describes the variables supplied in the dataset
metadata. Descriptions, units, and standard names are transcribed from
the provided metadata.

## Site, sampling, and geographic information

  --------------------------------------------------------------------------------------
  Variable                       Description       Units             Standard name
  ------------------------------ ----------------- ----------------- -------------------
  `Site_ID`                      Unique identifier Unitless          `site`
                                 for each site.                      

  `Sample_ID`                    Unique identifier Unitless          `sample`
                                 for each sampling                   
                                 event.                              

  `Data_Source`                  Source of         Unitless          `sampling_method`
                                 dataset.                            

  `Latitude_Degrees`             Latitude          Degrees North     `latitude`
                                 coordinates;                        
                                 positive values                     
                                 indicate North                      
                                 and negative                        
                                 values indicate                     
                                 South.                              

  `Longitude_Degrees`            Longitude         Degrees East      `longitude`
                                 coordinates;                        
                                 positive values                     
                                 indicate East and                   
                                 negative values                     
                                 indicate West.                      

  `Ocean_Name`                   Ocean in which    Unitless          `region`
                                 sampling took                       
                                 place.                              

  `Reef_ID`                      Unique identifier Unitless          `Site_ID`
                                 from Reef Check                     
                                 data.                               

  `Realm_Name`                   Realm as defined  Unitless          `region`
                                 by the Marine                       
                                 Ecoregions of the                   
                                 World (MEOW),                       
                                 Spalding et                         
                                 al. (2007).                         

  `Ecoregion_Name`               Ecoregion (150    Unitless          `region`
                                 ecoregions) as                      
                                 defined by Veron                    
                                 et al.                              

  `Country_Name`                 Country where     Unitless          `region`
                                 sampling took                       
                                 place.                              

  `State_Island_Province_Name`   State, territory  Unitless          `region`
                                 (e.g., Guam), or                    
                                 island group                        
                                 (e.g., Hawaiian                     
                                 Islands) where                      
                                 sampling took                       
                                 place.                              

  `City_Town_Name`               Region, city, or  Unitless          `region`
                                 nearest town                        
                                 where sampling                      
                                 took place.                         

  `Site_Name`                    Accepted name of  Unitless          `site`
                                 the site or name                    
                                 given by the team                   
                                 that sampled the                    
                                 reef.                               

  `Distance_to_Shore`            Distance of the   Meters (m)        `length`
                                 sampling site                       
                                 from the nearest                    
                                 land.                               

  `Exposure`                     Site exposure to  Unitless          `site_descrip`
                                 fetch. A site was                   
                                 considered                          
                                 exposed if it had                   
                                 more than 20 km                     
                                 of fetch, strong                    
                                 seasonal winds,                     
                                 or faced the                        
                                 prevailing winds.                   
                                 Otherwise, it was                   
                                 considered                          
                                 sheltered or                        
                                 "sometimes."                        
                                 "Sometimes"                         
                                 refers to sites                     
                                 with more than 20                   
                                 km of fetch                         
                                 through a narrow                    
                                 geographic window                   
                                 and therefore                       
                                 potentially                         
                                 exposed during                      
                                 cyclone seasons.                    

  `Cyclone_Frequency`            Number of cyclone Unitless          `site_descrip`
                                 events from 1964                    
                                 to 2014.                            

  `Depth_m`                      Depth of the      Meters (m)        `depth_w`
                                 sampling site.                      

  `Substrate_Name`               Type of substrate Unitless          `site_descrip`
                                 from Reef Check                     
                                 data.                               
  --------------------------------------------------------------------------------------

## Sampling date and comments

  ----------------------------------------------------------------------------
  Variable               Description       Units             Standard name
  ---------------------- ----------------- ----------------- -----------------
  `Date_Day`             Day of the        Unitless          `day`
                         sampling event.                     

  `Date_Month`           Month of the      Unitless          `month`
                         sampling event.                     

  `Date_Year`            Year of the       Unitless          `year`
                         sampling event.                     

  `Date`                 Date of sampling  Unitless          `date`
                         event in                            
                         `YYYY-MM-DD`                        
                         format.                             

  `Site_Comments`        Comments about    Unitless          `comment`
                         issues with the                     
                         site or                             
                         additional                          
                         information.                        

  `Sample_Comments`      Comments about    Unitless          `comment`
                         issues or                           
                         additional                          
                         information                         
                         concerning the                      
                         sampling event.                     

  `Bleaching_Comments`   Comments about    Unitless          `comment`
                         issues or                           
                         additional                          
                         information                         
                         concerning the                      
                         bleaching value.                    
  ----------------------------------------------------------------------------

## Reef condition and bleaching measurements

  ----------------------------------------------------------------------------
  Variable              Description       Units             Standard name
  --------------------- ----------------- ----------------- ------------------
  `Percent_Cover`       Average cover     Percent (%)       `percent_cover`
                        value.                              

  `Bleaching_Level`     Reef Check data   Unitless          `site_descrip`
                        indicating                          
                        whether the                         
                        observation                         
                        refers to a coral                   
                        population or                       
                        coral colony.                       

  `Percent_Bleaching`   Average of four   Percent (%)       `bleach_percent`
                        transect segments                   
                        (Reef Check) or                     
                        average of a                        
                        bleaching code.                     
  ----------------------------------------------------------------------------

## Environmental and temperature measurements

  -----------------------------------------------------------------------------------------------
  Variable                                  Description       Units             Standard name
  ----------------------------------------- ----------------- ----------------- -----------------
  `Turbidity`                               Kd490 with a      Reciprocal meters `turbidity`
                                            100-km buffer.    (m⁻¹)             
                                            Turbidity is                        
                                            positively                          
                                            related to the                      
                                            diffuse                             
                                            attenuation                         
                                            coefficient of                      
                                            light at the 490                    
                                            nm wavelength                       
                                            (Kd490), or the                     
                                            rate at which                       
                                            light at 490 nm                     
                                            is attenuated                       
                                            with depth. A                       
                                            Kd490 value of                      
                                            0.1 m⁻¹ means                       
                                            that light                          
                                            intensity is                        
                                            reduced by one                      
                                            natural-log value                   
                                            within 10 m of                      
                                            water. Higher                       
                                            Kd490 values                        
                                            represent higher                    
                                            attenuation and                     
                                            higher turbidity.                   

  `ClimSST`                                 Climatological    Degrees Celsius   `temp_ss`
                                            sea surface       (°C)              
                                            temperature (SST)                   
                                            based on weekly                     
                                            SSTs for the                        
                                            study time frame,                   
                                            created using a                     
                                            harmonics                           
                                            approach.                           

  `Temperature_Kelvin`                      Temperature in    Kelvin (K)        `temperature`
                                            Kelvin.                             

  `Temperature_Mean`                        Mean temperature. Degrees Celsius   `temperature`
                                                              (°C)              

  `Temperature_Minimum`                     Minimum           Degrees Celsius   `temperature`
                                            temperature.      (°C)              

  `Temperature_Maximum`                     Maximum           Degrees Celsius   `temperature`
                                            temperature.      (°C)              

  `Temperature_Kelvin_Standard_Deviation`   Standard          Kelvin (K)        `temperature`
                                            deviation of                        
                                            temperature.                        

  `Windspeed`                               Wind speed.       Meters per hour   `wind_speed`
                                                              (as supplied in   
                                                              metadata)         

  `SSTA`                                    Sea Surface       Degrees Celsius   `temp_ss`
                                            Temperature       (°C)              
                                            Anomaly: weekly                     
                                            SST minus weekly                    
                                            climatological                      
                                            SST.                                

  `SSTA_Standard_Deviation`                 Standard          Degrees Celsius   `temp_ss`
                                            deviation of      (°C)              
                                            weekly SST                          
                                            anomalies over                      
                                            the entire time                     
                                            period.                             

  `SSTA_Mean`                               Mean SSTA over    Degrees Celsius   `temp_ss`
                                            the entire time   (°C)              
                                            period.                             

  `SSTA_Minimum`                            Minimum SSTA over Degrees Celsius   `temp_ss`
                                            the entire time   (°C)              
                                            period.                             

  `SSTA_Maximum`                            Maximum SSTA over Degrees Celsius   `temp_ss`
                                            the entire time   (°C)              
                                            period.                             

  `SSTA_Frequency`                          Number of times   SSTA per time     `temp_ss`
                                            over the previous period            
                                            52 weeks that                       
                                            SSTA was greater                    
                                            than or equal to                    
                                            1 °C.                               

  `SSTA_Frequency_Standard_Deviation`       Standard          SSTA per time     `temp_ss`
                                            deviation of      period            
                                            SSTA_Frequency                      
                                            over the entire                     
                                            time period.                        

  `SSTA_FrequencyMax`                       Maximum           SSTA per time     `temp_ss`
                                            SSTA_Frequency    period            
                                            over the entire                     
                                            time period.                        

  `SSTA_FrequencyMean`                      Mean              SSTA per time     `temp_ss`
                                            SSTA_Frequency    period            
                                            over the entire                     
                                            time period.                        

  `SSTA_DHW`                                Sea Surface       Weeks             `time_elapsed`
                                            Temperature                         
                                            Degree Heating                      
                                            Weeks: sum of the                   
                                            previous 12 weeks                   
                                            when SSTA was                       
                                            greater than or                     
                                            equal to 1 °C.                      

  `SSTA_DHW_Standard_Deviation`             Standard          Weeks             `temp_ss`
                                            deviation of                        
                                            SSTA_DHW over the                   
                                            entire time                         
                                            period.                             

  `SSTA_DHWMax`                             Maximum SSTA_DHW  Weeks             `temp_ss`
                                            over the entire                     
                                            time period.                        

  `SSTA_DHWMean`                            Mean SSTA_DHW     Weeks             `temp_ss`
                                            over the entire                     
                                            time period.                        

  `TSA`                                     Thermal Stress    Degrees Celsius   `temperature`
                                            Anomaly: weekly   (°C)              
                                            sea surface                         
                                            temperature minus                   
                                            the maximum of                      
                                            weekly                              
                                            climatological                      
                                            sea surface                         
                                            temperature.                        

  `TSA_Standard_Deviation`                  Standard          Degrees Celsius   `temperature`
                                            deviation of TSA  (°C)              
                                            over the entire                     
                                            time period.                        

  `TSA_Minimum`                             Minimum TSA.      Degrees Celsius   `temperature`
                                                              (°C)              

  `TSA_Maximum`                             Maximum TSA.      Degrees Celsius   `temperature`
                                                              (°C)              

  `TSA_Mean`                                Mean TSA over the Degrees Celsius   `temperature`
                                            entire time       (°C)              
                                            period.                             

  `TSA_Frequency`                           Number of times   TSA per time      `temperature`
                                            over the previous period            
                                            52 weeks that TSA                   
                                            was greater than                    
                                            or equal to 1 °C.                   

  `TSA_Frequency_Standard_Deviation`        Standard          TSA per time      `temperature`
                                            deviation of      period            
                                            frequency of                        
                                            thermal stress                      
                                            anomalies over                      
                                            the entire time                     
                                            period.                             

  `TSA_FrequencyMax`                        Maximum           TSA per time      `temperature`
                                            TSA_Frequency     period            
                                            over the entire                     
                                            time period.                        

  `TSA_FrequencyMean`                       Mean              TSA per time      `temperature`
                                            TSA_Frequency     period            
                                            over the entire                     
                                            time period.                        

  `TSA_DHW`                                 Thermal Stress    Weeks             `time_elapsed`
                                            Anomaly (TSA)                       
                                            Degree Heating                      
                                            Weeks: sum of the                   
                                            previous 12 weeks                   
                                            when TSA was                        
                                            greater than or                     
                                            equal to 1 °C.                      

  `TSA_DHW_Standard_Deviation`              Standard          Weeks             `time_elapsed`
                                            deviation of                        
                                            TSA_DHW over the                    
                                            entire time                         
                                            period.                             

  `TSA_DHWMax`                              Maximum TSA_DHW   Weeks             `time_elapsed`
                                            over the entire                     
                                            time period.                        

  `TSA_DHWMean`                             Mean TSA_DHW over Weeks             `time_elapsed`
                                            the entire time                     
                                            period.                             
  -----------------------------------------------------------------------------------------------

## Notes for this analysis

-   Identifier fields such as `Site_ID`, `Sample_ID`, and `Reef_ID`
    identify sites or sampling records; they are not continuous
    environmental measurements.
-   Date fields are present in the source metadata. For an analysis that
    excludes the time-series aspect, leave `Date`, `Date_Day`,
    `Date_Month`, and `Date_Year` out of the analytical feature set.
-   `Percent_Bleaching` is the bleaching outcome. When using PCA to
    summarize environmental conditions, exclude it from the PCA input
    and use it afterward to compare bleaching groups.
-   Categorical descriptors such as `Ocean_Name`, `Exposure`,
    `Substrate_Name`, and `Bleaching_Level` require appropriate
    categorical analysis or encoding; they should not be treated as
    continuous numerical measurements.
-   Several SSTA and TSA fields summarize or count conditions over
    specified time windows. Excluding date columns does not remove the
    time-window information already represented in these derived
    variables.
-   The metadata lists `Windspeed` in meters per hour. Verify this unit
    against the original dataset documentation before interpreting it.
