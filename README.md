# motorwayData
This is the repository for the data described in Metz-Peeters (2022). It contains all the code required to construct the data set from the individual data sources. Furthermore, `/open_data/` contains open versions of the constructed data sets, in accordance with the individual data licenses (Status: January 31, 2023). These new datasets are provided under the Creative Commons license [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). 
`data.csv` contains the main analysis data set, `me_data_imputed.csv` contains the imputed maximum extent data set with some missing crash data, and `geo.shp` contains the network shape with IDs, so the dataset can be merged to it. 

A guide on how to construct the final dataset can be found in the Documentation.pdf 


## Data sources used
<table class="tg">
<thead>
  <tr>
    <th class="tg-p1nr">Provider</th>
    <th class="tg-p1nr">Data</th>
    <th class="tg-0pky">Variable(s)</th>
    <th class="tg-0pky">Links</th>
    <th class="tg-0pky">License</th>
    <th class="tg-0pky">Consequence</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">OpenStreetMap contributors</td>
    <td class="tg-0pky">OpenStreetMap</td>
    <td class="tg-sg5v">maxspeed_100, maxspeed_120, maxspeed_130, tunnel, bridge, no_shoulder, main_Entry, sec_Entry, main_Exit, sec_Exit, ramps, node_area, straight, right_turn, left_turn, total_turn, n_lanes, (ms_cond, ms_change, overtaking_ht)</td>
    <td class="tg-0pky"><a href="https://overpass-turbo.eu/" target="_blank" rel="noopener noreferrer"><span style="color:#905">Download interface for motorway network</span></a> and <a href="https://download.geofabrik.de/europe/germany.html" target="_blank" rel="noopener noreferrer"><span style="color:#905">downloads for side roads</span></a></td>
    <td class="tg-0pky"><a href="https://wiki.osmfoundation.org/wiki/Licence/Licence_and_Legal_FAQ" target="_blank" rel="noopener noreferrer"><span style="color:#905">License defined</span></a>, <a href="https://opendatacommons.org/licenses/odbl/1-0/" target="_blank" rel="noopener noreferrer"><span style="color:#905">Open Data Commons Open Database License (ODbL) v1.0</span></a></td>
    <td class="tg-0pky">Open use under attribution</td>
  </tr>
  <tr>
    <td class="tg-0pky">Statistische Ämter des Bundes und der Länder</td>
    <td class="tg-0pky">Unfallatlas (locations of injury crashes)</td>
    <td class="tg-sg5v">total, fatal, severly_injured, lighlty_injured</td>
    <td class="tg-0pky"><a href="https://unfallatlas.statistikportal.de/_opendata2022.html" target="_blank" rel="noopener noreferrer"><span style="color:#905">Downloads</span></a></td>
    <td class="tg-0pky"><a href="https://unfallatlas.statistikportal.de/_opendata2022.html" target="_blank" rel="noopener noreferrer"><span style="color:#905">License defined</span></a>, <a href="https://www.govdata.de/dl-de/by-2-0" target="_blank" rel="noopener noreferrer"><span style="color:#905">Datenlizenz Deutschland – Namensnennung – Version 2.0</span></a></td>
    <td class="tg-0pky">Open use under attribution</td>
  </tr>
  <tr>
    <td class="tg-0pky">Bundesanstalt für Straßenwesen (BASt, Federal Highway Research Institute)</td>
    <td class="tg-0pky">Bundesfernstraßennetz (BISStra) (Federal highway network)</td>
    <td class="tg-sg5v">location</td>
    <td class="tg-0pky"><a href="https://www.bast.de/DE/Verkehrstechnik/Fachthemen/Daten/Daten-BISStra.html?nn=1817946" target="_blank" rel="noopener noreferrer"><span style="color:#905">Download</span></a></td>
    <td class="tg-0pky"><a href="https://www.govdata.de/web/guest/suchen/-/details/datensatz-bundesfernstrassennetzcc7e0" target="_blank" rel="noopener noreferrer"><span style="color:#905">License defined</span></a>, <a href="http://www.gesetze-im-internet.de/geonutzv/" target="_blank" rel="noopener noreferrer"><span style="color:#905">GeoNutzV</span></a></td>
    <td class="tg-0pky">Open use under attribution</td>
  </tr>
  <tr>
    <td class="tg-0pky">BASt</td>
    <td class="tg-0pky">Automatische Zählstellen 2017-2019 (automatic counting stations)</td>
    <td class="tg-sg5v" rowspan="2"><span style="font-weight:normal">AADT, AADT_HT, HT_share, MSV50, HT_share_MSV, AADT_day, AADT_night, night_share, sunday, holiday</span><br></td>
    <td class="tg-0pky"><a href="https://www.bast.de/DE/Verkehrstechnik/Fachthemen/v2-verkehrszaehlung/zaehl_node.html;jsessionid=9BE1F5EC97952EEFAB21C06D3F0BFD7B.live11314" target="_blank" rel="noopener noreferrer"><span style="color:#905">Downloads</span></a></td>
    <td class="tg-0pky"><a href="https://www.mcloud.de/web/guest/suche/-/results/detail/6CD31C11-50E5-4DB3-A7C7-8CA9774B525B" target="_blank" rel="noopener noreferrer"><span style="color:#905">License defined</span></a>, <a href="https://www.govdata.de/dl-de/by-nc-1-0" target="_blank" rel="noopener noreferrer"><span style="color:#905">Datenlizenz Deutschland Namensnennung - Version 1.0?</span></a></td>
    <td class="tg-0pky">Open use under attribution</td>
  </tr>
  <tr>
    <td class="tg-0pky">BASt</td>
    <td class="tg-0pky">Fortschreibung der SVZ 2015 und der termporären Messung 2016 bis 2019 auf das Jahr 2019 (Extrapolation of traffic counts of 2015 and of temporary measurements of 2016 to 2019.)</td>
    <td class="tg-0pky"><a href="https://www.bast.de/DE/Statistik/Verkehrsdaten/Manuelle-Zaehlung.html" target="_blank" rel="noopener noreferrer"><span style="color:#905">Download</span></a></td>
    <td class="tg-0pky">Unclear, waiting for response</td>
    <td class="tg-0pky">Unclear</td>
  </tr>
  <tr>
    <td class="tg-0pky">NASA</td>
    <td class="tg-0pky">Shuttle Radar Topography Mission (SRTM) 1 arc secon</td>
    <td class="tg-sg5v">down_change, up_change, max_slope, mean_slope, elevation</td>
    <td class="tg-0pky"><a href="https://www2.jpl.nasa.gov/srtm/" target="_blank" rel="noopener noreferrer"><span style="color:#905">Information</span></a></td>
    <td class="tg-0pky"><a href="https://www.earthdata.nasa.gov/learn/use-data/data-use-policy" target="_blank" rel="noopener noreferrer"><span style="color:#905">License defined and text</span></a></td>
    <td class="tg-0pky">Open use</td>
  </tr>
  <tr>
    <td class="tg-0pky">Deutscher Wetterdienst   (DWD, German Weather Service)</td>
    <td class="tg-0pky">Climate Data Center, annual grids</td>
    <td class="tg-sg5v">air_temp, frost_days, ice_days, snowcov_days, precip10mm, precip20mm, precip30mm, precipitation, summer_days, sunshine_dur, wind</td>
    <td class="tg-0pky"><a href="https://opendata.dwd.de/climate_environment/CDC/grids_germany/annual/" target="_blank" rel="noopener noreferrer"><span style="color:#905">Downloads</span></a></td>
    <td class="tg-0pky"><a href="https://opendata.dwd.de/climate_environment/CDC/Terms_of_use.pdf" target="_blank" rel="noopener noreferrer"><span style="color:#905">License defined</span></a>, <a href="http://www.gesetze-im-internet.de/geonutzv/" target="_blank" rel="noopener noreferrer"><span style="color:#905">GeoNutzV</span></a></td>
    <td class="tg-0pky">Open use under attribution</td>
  </tr>
  <tr>
    <td class="tg-0pky">Bundesamt für Bauwesen und Raumordnung (BBSR, Federal Institute for Research on Building; Urban Affairs and Spatial Development)</td>
    <td class="tg-0pky">INKAR - Indikatoren und Karten zur Raum- und Stadtentwicklung (indicators and maps of spatial and urban development)</td>
    <td class="tg-sg5v" rowspan="2">pop_dens, emp_quo, pop_18_25, pop_ol_65<br>fem_share, hh_inc, pcar_dens, gdp_p_cap, rurality</td>
    <td class="tg-0pky"><a href="https://www.inkar.de/" target="_blank" rel="noopener noreferrer"><span style="color:#905">Downloads</span></a></td>
    <td class="tg-0pky"><a href="https://www.bbsr.bund.de/BBSR/DE/service/nutzungshinweise/_node.html" target="_blank" rel="noopener noreferrer"><span style="color:#905">License defined</span></a>, <a href="https://www.govdata.de/dl-de/by-2-0" target="_blank" rel="noopener noreferrer"><span style="color:#905">Datenlizenz Deutschland – Namensnennung – Version 2.0</span></a></td>
    <td class="tg-0pky">Open use under attribution</td>
  </tr>
  <tr>
    <td class="tg-0pky">Bundesamt für Kartographie und Geodäsie (BKG, Federal Agency for Cartography and Geodesy)</td>
    <td class="tg-0pky">Verwaltungsgebiete 1:250 000 mit Einwohnerzahlen (Ebenen), Stand 31.12. (VG250-EW 31.12.)  (Administrative areas 1:250 000 with population numbers)</td>
    <td class="tg-0pky"><a href="https://gdz.bkg.bund.de/index.php/default/open-data/verwaltungsgebiete-1-250-000-mit-einwohnerzahlen-stand-31-12-vg250-ew-31-12.html" target="_blank" rel="noopener noreferrer"><span style="color:#905">Download</span></a></td>
    <td class="tg-0pky"><a href="https://gdz.bkg.bund.de/index.php/default/open-data/verwaltungsgebiete-1-250-000-mit-einwohnerzahlen-stand-31-12-vg250-ew-31-12.html" target="_blank" rel="noopener noreferrer"><span style="color:#905">License defined in tab "Nutzungsbedingungen"</span></a>, <a href="https://www.govdata.de/dl-de/by-2-0" target="_blank" rel="noopener noreferrer"><span style="color:#905">Datenlizenz Deutschland – Namensnennung – Version 2.0</span></a></td>
    <td class="tg-0pky">Open use under attribution</td>
  </tr>
  <tr>
    <td class="tg-0pky">BASt</td>
    <td class="tg-0pky">Zustandserfassung und -bewertung (ZEB) (Road condition measures)</td>
    <td class="tg-sg5v">asphalt, sub_mean, sub_var, sub_max, perf_mean, perf_var, per_max</td>
    <td class="tg-0pky">Received upon request. <a href="https://bmdv.bund.de/SharedDocs/DE/Artikel/StB/zustandserfassung-und-bewertung.html" target="_blank" rel="noopener noreferrer"><span style="color:#905">Information</span></a> only available in German</td>
    <td class="tg-0pky">Individual data usage agreement prohibits sharing the data.</td>
    <td class="tg-0pky">Excluded from new data set. Code to merge this data to the open data set will be provided.</td>
  </tr>
</tbody>
</table>



## References
Metz-Peeters, Maike (2022) : The effects of mandatory speed limits on crash frequency: A causal machine learning approach, Ruhr Economic Papers, 982.
