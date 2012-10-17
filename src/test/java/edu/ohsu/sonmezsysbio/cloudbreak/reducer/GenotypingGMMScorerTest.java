package edu.ohsu.sonmezsysbio.cloudbreak.reducer;

import edu.ohsu.sonmezsysbio.cloudbreak.ReadGroupInfo;
import edu.ohsu.sonmezsysbio.cloudbreak.io.ReadPairInfo;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: cwhelan
 * Date: 10/12/12
 * Time: 10:11 AM
 */

/**
 * These tests are mostly direct ports from my R prototype
 */
public class GenotypingGMMScorerTest {

    @Test
    public void testLogsumexp() throws Exception {
        double[] x = new double[] {-7.018195, -309.17497};
        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        assertEquals(-7.018195, scorer.logsumexp(x), 0.00001);
    }

    @Test
    public void testNnclean() throws Exception {
        double[] y = new double[] {260.0736, 197.4272,   194.8618,  1217.8588,
                1228.2190,  1151.7017,  4511.5326, 19719.9700, 19707.2091, 16788.1891,
                22556.1203,  9909.8687, 13709.7259,  9219.4829, 13076.1122,  8140.4713};
        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        double[] nonnoise = new double[] {260.0736, 197.4272,   194.8618,  1217.8588,
                1228.2190,  1151.7017};
        assertArrayEquals(nonnoise, scorer.nnclean(y, 30, 2), 0.000001);
    }

    @Test
    public void testNncleanNotEnoughValues() throws Exception {
        double[] y = new double[] {260.0736, 1217.8588};
        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        double[] nonnoise = new double[] {};
        assertArrayEquals(nonnoise, scorer.nnclean(y, 30, 2), 0.000001);
    }

    @Test
    public void testEstimateW() throws Exception {
        double[] y = new double[] {260.0736, 197.4272,   194.8618,  1217.8588,
                1228.2190,  1151.7017,  4511.5326, 19719.9700, 19707.2091, 16788.1891,
                22556.1203,  9909.8687, 13709.7259,  9219.4829, 13076.1122,  8140.4713};

        double sigma = 30;
        double[] initialW = new double[] {Math.log(.5),Math.log(.5)};

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        assertEquals(.5, scorer.estimateW(y, initialW, 200, sigma), 0.0001);
    }

    public void testEstimateWNoValues() throws Exception {
        double[] y = new double[] {};

        double sigma = 30;
        double[] initialW = new double[] {Math.log(.5),Math.log(.5)};

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        assertEquals(1, scorer.estimateW(y, initialW, 200, sigma), 0.0001);
    }

    @Test
    public void testReduce_2_91700() throws Exception {
        double[] y = new double[] {152,
                216,
                194,
                169,
                202,
                237,
                178,
                202,
                208,
                210,
                247,
                227,
                182,
                207,
                191,
                147,
                251,
                236,
                209,
                200,
                204,
                184,
                158,
                152,
                255,
                227,
                185,
                183,
                204,
                226,
                276,
                207,
                223
        };

        ReadGroupInfo readGroupInfo = new ReadGroupInfo();
        readGroupInfo.isize = 200;
        readGroupInfo.isizeSD = 30;

        List<ReadPairInfo> rpis = new ArrayList<ReadPairInfo>();
        for (int i = 0; i < y.length; i++) {
            double isize = y[i];
            ReadPairInfo rpi = new ReadPairInfo();
            rpi.insertSize = (int) isize;
            rpi.readGroupId = 0;
            rpi.pMappingCorrect = 0;
            rpis.add(rpi);
        }

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();

        Map<Short, ReadGroupInfo>  rgis = new HashMap<Short, ReadGroupInfo>();
        rgis.put((short) 0, readGroupInfo);
        double score = scorer.reduceReadPairInfos(rpis.iterator(), rgis);
        assertEquals(1, score, 0.00001);
    }

    @Test
    public void testReduce_homdel_100() throws Exception {
        double[] y = new double[] {152,
                216,
                194,
                169,
                202,
                237,
                178,
                202,
                208,
                210,
                247,
                227,
                182,
                207,
                191,
                147,
                251,
                236,
                209,
                200,
                204,
                184,
                158,
                152,
                255,
                227,
                185,
                183,
                204,
                226,
                276,
                207,
                223
        };

        ReadGroupInfo readGroupInfo = new ReadGroupInfo();
        readGroupInfo.isize = 100;
        readGroupInfo.isizeSD = 15;

        List<ReadPairInfo> rpis = new ArrayList<ReadPairInfo>();
        for (int i = 0; i < y.length; i++) {
            double isize = y[i];
            ReadPairInfo rpi = new ReadPairInfo();
            rpi.insertSize = (int) isize;
            rpi.readGroupId = 0;
            rpis.add(rpi);
        }

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();

        Map<Short, ReadGroupInfo>  rgis = new HashMap<Short, ReadGroupInfo>();
        rgis.put((short) 0, readGroupInfo);
        double score = scorer.reduceReadPairInfos(rpis.iterator(), rgis);
        assertEquals(0, score, 0.0001);
    }

    @Test
    public void testLikelihood() throws Exception {
        double[] y = new double[] {152,
                216,
                194,
                169,
                202,
                237,
                178,
                202,
                208,
                210,
                247,
                227,
                182,
                207,
                191,
                147,
                251,
                236,
                209,
                200,
                204,
                184,
                158,
                152,
                255,
                227,
                185,
                183,
                204,
                226,
                276,
                207,
                223
        };

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();
        double l = scorer.likelihood(y, new double[] { Math.log(.5), Math.log(.5)}, new double[] {100, 1100}, 15);
        assertEquals(-908.0622, l, 0.0001);
    }

    @Test
    public void testReduce_2_16168700() throws Exception {
        double[] y = new double[] {492,549,558,517,562,523,533,537,544,546,540,506,516,519,560,515,498,546,493,
                518,520,517,153,164,189,208,154,203,540,194,497,566,551,555,139,491,510,229,199,190,213,190,151,
                210,218,226,186,173,171,226,212,181,216,146,171,184,214,199,192,184,136,208,214,226,203,200,217,
                176,174,222,182,186,188,202,211,190,216,163,201,244,241,216,254,189,189,158,213,240,203,238,241,
                178,195,196,213,225,187,239,219,171,235,227,216,163,214,220,195,180,189,206,211,178,206,207,201,
                205,203,197,221,215,215,161,199,136,203,212,172,211,222,180,157,207,187,202,188,191,180,197,197,
                199,209,237,255,174,187,212,215,218,142,201,150,155,206,192,182,156,201,248,195,206,165,178,217,
                184,178,164,167,130,196,195,139,211,191,176,212,195,149,226,199,247,168,223,176,198,196,226,144,
                196,213,191,203,173,188,186,265,222,219,209,218,206,226,208,188,240,172,210,198,155,206,163,183,
                165,197,168,135,225,226,190,222,189,152,245,212,231,193,150,162,210,185,217,244,233,167,100,169,
                187,203,141,189,214,158,154,164,210,220,176,208,184,218,157,171,184,203,200,189,187,223,199,174,
                194,210,155,145,189,250,209,174,196,165,262,219,182,196,206,201,180,168,212,187,215,207,243,170,
                195,210,146,244,271,134,147,220,206,206,214,167,192,173,181,185,234,189,207,218,185,221,165,244,
                247,251,195,170,204,183,183,201,202,169,183,243,186,202,207,195,134,146,165,214,217,180,187,166,
                219,179,204,190,172,204,159,247,217,234,215,196,216,178,235,190,225,181,228,243,201,215,239,231,
                232,215,150,224,189,212,192,208,237,189,155,279,182,214,200,176,211,236,231,208,177,235,207,221,
                189,222,209,241,191,226,232,217,179,208,193,217,164,218,192,186,147,190,255,203,187,197,159,179,
                217,214,191,187,210,194,205,161,188,147,204,165,170,191,170,254,170,187,149,206,209,180,172,277,
                219,245,198,206,181,115,235,188,146,219,251,176,196,221,252,238,217,164,220,213,183,200,169,140,
                182,201,259,175,179,245,173,208,184,171,210,161,185,197,225,229,162,228,239,151,176,175,191,185,
                172,210,222,196,244,155,210,226,206,202,226,182,203,211,54,183,162,218,178,222,186,175,185,206
        };

        double[] q = new double[]{-2.1696611175533603E-10, -6.785146088283706E-8, -8.024704012413976E-8, -1.7379881135464118E-4, -0.6931471805599502, -0.693147180559951, -0.6931471805678728, -0.693147180572517, -0.6931471805725841, -0.6931471805798182, -0.6931471805997381, -0.693147180622788, -0.6931471806495079, -0.6931471807594544, -0.6931471845248532, -0.6931475377340768, -0.6931503533452836, -0.6931512316685288, -0.6931840359492512, -0.693305956786865, -0.6937753785184062, -5.64897333105462, -7.084857412316782, -7.354808352597096, -7.359377330663584, -7.872364521913692, -8.04267742026059, -8.669566545094245, -8.969160504596012, -9.221512806273967, -9.443121607473534, -9.662320273839363, -9.89927891367327, -9.899475907902985, -10.088747376089458, -10.136296317056333, -10.390994468392968, -10.466896786383858, -10.54685718548888, -10.563884362458577, -10.581221493505096, -10.62770119673156, -10.63334549133029, -10.646393514132974, -10.657253069096956, -10.663702915211145, -10.676735018643225, -10.681529869083171, -10.686913818512526, -10.699011191869902, -10.712076390390033, -10.71494371186775, -10.7264211634449, -10.740214964542467, -10.749791882491355, -10.750311639419573, -10.754497756231224, -10.756823839754684, -10.761609867517738, -10.766403150113543, -10.771343223309085, -10.775062291368782, -10.77835101998959, -10.78197573924654, -10.790276905621603, -10.795256866828304, -10.80181701129367, -10.803087896128151, -10.80314009117173, -10.80380175386974, -10.814394966508459, -10.81537122505761, -10.817492173087448, -10.827675844546354, -10.829443481011559, -10.830179309986828, -10.838212544645792, -10.844131294622137, -10.848542068923848, -10.853854857544324, -10.856020915168198, -10.859558263487916, -10.860158700196953, -10.880641680567678, -10.888588199397667, -10.88874066031219, -10.91229133295613, -10.924107483771525, -10.928346113518579, -10.92964687162989, -10.932350396150245, -10.947241656379116, -10.960930793781891, -10.968987987336805, -10.97537609031065, -10.97888333176698, -10.983714693620161, -10.98396459725721, -10.988026161635448, -11.000349775425226, -11.005190098864151, -11.008058181344321, -11.019278024115216, -11.021291910823129, -11.028322692242273, -11.029149326715867, -11.061707802901795, -11.068395939657464, -11.069852995949653, -11.069933889625611, -11.076805056632331, -11.083382084815497, -11.08509571105233, -11.089824851323177, -11.089885487887706, -11.094113103444364, -11.100555423526552, -11.10133122503366, -11.113614986629823, -11.120602422354539, -11.145674693493445, -11.155171514432995, -11.155584476787158, -11.16235322235817, -11.1743979194754, -11.179930631664496, -11.18119861073529, -11.186147109040785, -11.186253203506407, -11.194110915846878, -11.199850352222427, -11.205512744495765, -11.206608319859358, -11.208810637615429, -11.208935886080297, -11.212548365317536, -11.21424709756695, -11.23214588034072, -11.23570216013649, -11.241624346768159, -11.244808957153834, -11.266084963986366, -11.268474719034423, -11.274459210338485, -11.277656936264519, -11.287998772134518, -11.29176173512852, -11.32731360798145, -11.3311640427555, -11.333039066255113, -11.337135957453896, -11.358383144232379, -11.366000617155429, -11.36659303658937, -11.384491692834233, -11.393971852444107, -11.397996304813965, -11.398666888576694, -11.402682214394648, -11.405528284621504, -11.407568257147572, -11.407675455579561, -11.41065714680622, -11.413605020736467, -11.4148346590077, -11.415836742763219, -11.420320556150081, -11.428985678366118, -11.433031937419411, -11.434169066321804, -11.435919487005071, -11.440941652690302, -11.448930587418992, -11.452612464056163, -11.453244486778988, -11.456067248558421, -11.459937035668403, -11.465816007733807, -11.4667715838359, -11.469548419058178, -11.506416283676629, -11.507980892477075, -11.508295681284594, -11.510764489174388, -11.521049706794319, -11.534339191667677, -11.554287990779184, -11.554403010730585, -11.557098234728846, -11.558454457201158, -11.558527724793233, -11.56231529416586, -11.571770761162302, -11.571770950887338, -11.607173819094001, -11.607813409351401, -11.609557604640045, -11.612361592929172, -11.616746254410995, -11.617624768597299, -11.619735777410394, -11.627618299073958, -11.62944394218113, -11.632677032446942, -11.633634976353044, -11.656295440568524, -11.67754449036551, -11.679715107322874, -11.681644956169858, -11.683492771248496, -11.684179635522574, -11.688069451670032, -11.69307353975453, -11.705608623096106, -11.718261062653541, -11.721859752278295, -11.725105691811862, -11.74293467318785, -11.74604369023696, -11.750169688924437, -11.755383243031819, -11.75838466961116, -11.760354632643484, -11.765392247336596, -11.773565769921127, -11.795011574718444, -11.798918108029099, -11.80044509306538, -11.802451468952714, -11.802870247278562, -11.808387690341426, -11.811275898644656, -11.81515379404222, -11.823109810566377, -11.830430003357971, -11.831698794922358, -11.833383051107955, -11.845948060608585, -11.857525221947936, -11.872834546949225, -11.884868237981923, -11.892365054656624, -11.89518643359543, -11.906122690678668, -11.907464313404134, -11.919426054164733, -11.936547937965024, -11.948251494311359, -11.960643538424119, -11.964765517009852, -11.971871206040912, -11.979188840635903, -11.979872541286628, -11.980059958768567, -11.983324560267182, -11.989347998124508, -11.990776747102908, -11.99889265952341, -12.008243796387097, -12.023278749706506, -12.031611349208564, -12.033054852806274, -12.049962180295484, -12.051187384115412, -12.05373847617144, -12.059378026508043, -12.061165014529005, -12.08425921677133, -12.09290153858423, -12.096657564011192, -12.104708475256558, -12.105411877986661, -12.120906812725586, -12.125078259994385, -12.128819793798097, -12.131577706216952, -12.137083478593583, -12.139965138023623, -12.144840271962934, -12.16647882827111, -12.176385324600925, -12.19331777410186, -12.206280608871772, -12.206368172521795, -12.21221581648059, -12.213997443665919, -12.217147145466207, -12.223044019810153, -12.227008013842497, -12.240297236766862, -12.246202784275956, -12.247633359841192, -12.249949482076318, -12.252280350545725, -12.259958059662193, -12.269491325434306, -12.271041060855865, -12.271560124839013, -12.278694185360683, -12.27872132138998, -12.281563915928338, -12.28822897720492, -12.296341433158773, -12.308269982919533, -12.327691266877867, -12.336869339082376, -12.341056064100963, -12.356900031894881, -12.376800778971106, -12.380930491424309, -12.399361188991575, -12.400429215097255, -12.400505214268058, -12.415877378043504, -12.418190646488394, -12.428500964237001, -12.434414081290308, -12.454876625161518, -12.456919070770772, -12.469664441900012, -12.493521693362812, -12.494697985905036, -12.504294353427472, -12.520024626960016, -12.535297564869646, -12.539342819183414, -12.540348346241775, -12.554023456234265, -12.559622674013278, -12.571204825571195, -12.574741730551908, -12.583825658741782, -12.626304184802747, -12.627749295535345, -12.637190715097805, -12.647945005733536, -12.669042955168685, -12.675206250932629, -12.681832925909207, -12.687801255687505, -12.691598881155024, -12.711106689629208, -12.73904118528932, -12.742302306626886, -12.74955717594839, -12.765080380251657, -12.768093773024798, -12.775903088506945, -12.823522421893, -12.826894576679015, -12.827617882552838, -12.831061272166139, -12.831646753233883, -12.843083873122257, -12.86085924799955, -12.903402521217178, -12.906104025278378, -12.906447168159744, -12.910875138821364, -12.91455067260042, -12.934919274290205, -12.956163559665734, -12.963918515173502, -12.97051594259994, -12.9753163463167, -12.984005448337019, -12.999829630789407, -13.045536310822637, -13.048772683530792, -13.052781597159033, -13.085228846958847, -13.08691508819893, -13.093911862249723, -13.142853842381857, -13.152551785535437, -13.154477154375927, -13.161777998003828, -13.167025161781197, -13.169691464089375, -13.181291817433694, -13.182216311875594, -13.219039826262367, -13.219499341503145, -13.227825837956352, -13.232727488905265, -13.252450113577616, -13.25439216258476, -13.264044774348491, -13.287510294333067, -13.297580706956184, -13.323571082504627, -13.334592054896639, -13.360688679044983, -13.376950374813497, -13.385030581320638, -13.400183248891341, -13.445952135319157, -13.449928892501957, -13.451586983108424, -13.469723421966817, -13.522936682874917, -13.551028122190452, -13.557758633290474, -13.57476024509103, -13.605446257818688, -13.606780427143324, -13.629226177962948, -13.637966053406343, -13.66253085031366, -13.712552861438235, -13.741313735475032, -13.791349633886988, -13.800682412370243, -13.802428810557252, -13.825743843031923, -13.83659679813952, -13.844631092429726, -13.846398874836751, -13.8602601997209, -13.861075786933808, -13.873522177022677, -13.8830254730683, -13.954565735702303, -13.962626346683482, -14.011017965921496, -14.021602339057813, -14.061349097166467, -14.1346480370544, -14.18642244862005, -14.264723939869263, -14.27575260774345, -14.292091994883386, -14.317414739954774, -14.32918558019255, -14.3363631404843, -14.40035941836633, -14.435348487784918, -14.443018727672229, -14.472990208844504, -14.500717254745936, -14.610011834805242, -14.616852240862958, -14.666163986368062, -14.684036947144627, -14.73093485648133, -14.735052619883703, -14.760689339795253, -14.763644797399394, -14.82090441612769, -14.848015393628826, -14.949915735929597, -15.004929672272947, -15.029480148258049, -15.032346096384781, -15.075222638833846, -15.080903343346119, -15.097314551862558, -15.115440735308157, -15.131922627917003, -15.226072618784602, -15.22890069250707, -15.233054064019132, -15.269549783702754, -15.303170788599967, -15.303788808882992, -15.416808319958273, -15.56464735941875, -15.796329199467257, -15.846777242562403, -16.005857824831633, -16.12351236735087, -16.150061512413757, -16.247090127882622, -16.278741427053546, -16.31169763271048, -16.416623168840637, -16.532872414202522, -16.684488046806592, -16.713654174537385, -16.87719518865177, -17.008598037470506, -17.10448955459016, -17.123895973086224, -17.25176272912007, -17.29356856286769, -17.409667809368557, -17.928570585634525, -17.960983028162417, -18.418115592931766, -18.420680743952364, -18.420680743952364, -18.420680743952364, -18.420680743952364, -18.45900866840194, -18.607842667340645, -19.147320224015054, -19.35446227099996, -19.538798616221108, -19.600596948545835};

        ReadGroupInfo readGroupInfo = new ReadGroupInfo();
        readGroupInfo.isize = 200;
        readGroupInfo.isizeSD = 30;

        List<ReadPairInfo> rpis = new ArrayList<ReadPairInfo>();
        for (int i = 0; i < y.length; i++) {
            double isize = y[i];
            ReadPairInfo rpi = new ReadPairInfo();
            rpi.insertSize = (int) isize;
            rpi.readGroupId = 0;
            rpi.pMappingCorrect = q[i];
            rpis.add(rpi);
        }

        GenotypingGMMScorer scorer = new GenotypingGMMScorer();

        Map<Short, ReadGroupInfo>  rgis = new HashMap<Short, ReadGroupInfo>();
        rgis.put((short) 0, readGroupInfo);
        double score = scorer.reduceReadPairInfos(rpis.iterator(), rgis);
        assertEquals(0, score, 0.0001);
    }

}
