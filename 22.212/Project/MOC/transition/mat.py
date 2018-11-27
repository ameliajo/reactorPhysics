


class Material:
    def __init__(self,SigT,nuSigF,SigS_matrix,chi):
        self.SigT = SigT
        self.nuSigF = nuSigF
        self.SigS_matrix = SigS_matrix
        self.chi = chi




modTotal = [2.595327130073840638e-01, 7.958026873794312728e-01, 1.296887889543246963e+00, 1.301285995868290524e+00, 1.305642831365466217e+00, 1.351106574249372638e+00, 1.501539480826580730e+00, 1.753774326215949753e+00, 2.116169195925658997e+00, 3.055718902416071803e+00]

modChi = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
modNuFission = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


modScatter = [[1.600173951507630754e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],[9.891213788419001007e-02, 6.702401738837177048e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],[5.954104341416693305e-04, 1.239787164433226729e-01, 1.073196277535670040e+00, 2.176095785464819991e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],[6.727801515725074472e-06, 7.047290227499326560e-04, 1.674932738620249495e-01, 6.666034219614380696e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],[0.0, 1.366584796967549795e-04, 3.330685450284268101e-02, 3.854518464793834642e-01, 5.372684488027054117e-01, 3.308903403507236902e-05, 0.0, 0.0, 0.0, 0.0],[0.0, 6.430987279847292354e-05, 1.851083370803453462e-02, 2.094927412666981736e-01, 6.507637689781351487e-01, 8.513808457224117943e-01, 3.886172332861179395e-03, 0.0, 0.0, 0.0],[0.0, 2.679578033269704583e-06, 1.849502599351017431e-03, 2.093404145617158013e-02, 6.639350959489830062e-02, 2.984630869963527666e-01, 7.719221759520857873e-01, 2.233489487494325423e-02, 0.0, 0.0], [0.0, 0.0, 7.376933444704911799e-04, 9.009036551824358557e-03, 2.623527899356360021e-02, 1.137766435295963358e-01, 4.422749338453474954e-01, 9.362708164094512009e-01, 7.594513745234211799e-02, 9.613331371229879435e-03],[0.0, 2.679578033269704583e-06, 3.846543867596133293e-04, 4.917976475150495456e-03, 1.431665975092147498e-02, 5.669805981909651488e-02, 1.978382594040979381e-01, 5.824687986359826652e-01, 1.431878723565403355e+00, 3.041259809858234964e-01],[0.0, 2.679578033269704583e-06, 2.002310506419904801e-04, 2.807163563249618493e-03, 8.232079356779846771e-03, 2.508148779858486443e-02, 7.526339261165088113e-02, 1.964606601276629727e-01, 5.835118060921619110e-01, 2.699815261878973693e+00]]




fuelTotal = [2.743855399193625977e-01, 4.322740579381196824e-01, 5.515839976044424331e-01, 6.129862382792780062e-01, 6.357097773893524151e-01, 4.334822130047676914e-01, 5.147607488413211696e-01, 5.849240737377591204e-01, 6.736396988572106448e-01, 9.546312458085901564e-01]
fuelChi = [7.627512073710626117e-01, 2.370893283062106005e-01, 1.594643227269019819e-04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
fuelNuFission = [2.656842845690966792e-02, 2.992139538647185012e-03, 2.775219454267422239e-02, 7.818989319222159617e-02, 6.626235844805078434e-02, 6.676923941837201171e-02, 2.160541725000551527e-01, 3.352685286318000668e-01, 5.080873787995320301e-01, 1.026838032506582898e+00]

fuelScatter = [[2.144511894524181728e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [4.888532775969586797e-02, 4.229112680264441315e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [7.079699892787236853e-06, 2.545554182608533549e-03, 4.987570140671064300e-01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 4.845188514006277050e-03, 4.334909573128676863e-01, 9.150478812363896559e-05, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.781913501659770915e-02, 3.692218200788832116e-01, 1.157811517168803189e-04, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 2.836648431832809042e-02, 3.721592153352922261e-01, 2.860542455146882629e-03, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 1.316045757848540235e-02, 3.574836732917382465e-01, 7.986282511608915752e-03, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.079289819363996775e-02, 3.427446244565491296e-01, 1.206318284097213783e-02, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.682672032439341814e-04, 4.459007735648311960e-02, 3.430283730419118049e-01, 2.829193366077875985e-02], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.319044282925950294e-05, 3.898467625436117384e-02, 3.760159297174331949e-01]]




modClass  = Material(modTotal,modNuFission,modScatter,modChi)
fuelClass = Material(fuelTotal,fuelNuFission,fuelScatter,fuelChi)


