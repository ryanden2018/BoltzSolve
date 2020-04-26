var gaussQuad =[
    [ 0.1392518728556320, -0.0697392733197222],
    [ 0.1392518728556320, 0.0697392733197222],
    [ 0.1365414983460152, -0.2078604266882213],
    [ 0.1365414983460152, 0.2078604266882213],
    [ 0.1311735047870624, -0.3419358208920842],
    [ 0.1311735047870624, 0.3419358208920842],
    [ 0.1232523768105124, -0.4693558379867570],
    [ 0.1232523768105124,  0.4693558379867570],
    [ 0.1129322960805392, -0.5876404035069116],
    [ 0.1129322960805392, 0.5876404035069116],
    [ 0.1004141444428810, -0.6944872631866827],
    [ 0.1004141444428810, 0.6944872631866827],
    [ 0.0859416062170677, -0.7878168059792081],
    [ 0.0859416062170677, 0.7878168059792081],
    [ 0.0697964684245205, -0.8658125777203002],
    [ 0.0697964684245205, 0.8658125777203002],
    [ 0.0522933351526833, -0.9269567721871740],
    [ 0.0522933351526833, 0.9269567721871740],
    [ 0.0337749015848142, -0.9700604978354287],
    [ 0.0337749015848142, 0.9700604978354287],
    [ 0.0146279952982722, -0.9942945854823992],
    [ 0.0146279952982722, 0.9942945854823992]
];

function hermiteEval(x,n)
{
    var a = Math.exp(-0.5*x*x) * Math.pow(3.141592653589793238462,-0.25);
    var b = Math.sqrt(2.0) * x * Math.exp(-0.5*x*x) * Math.pow(3.141592653589793238462,-0.25);
    for(var i = 0; i < n; i++)
    {
        var c =  Math.sqrt(2.0) * Math.pow((i+2.0),-0.5) * x * b - Math.pow(i+1.0,0.5) * Math.pow(i+2.0,-0.5) * a;
        a = b;
        b = c;
    }
    return a;
}

function hermiteDerivEval(x,n)
{
    var a = Math.exp(-0.5*x*x) * Math.pow(3.141592653589793238462,-0.25);
    var b = Math.sqrt(2.0) * x * Math.exp(-0.5*x*x) * Math.pow(3.141592653589793238462,-0.25);
    var c = x * b - Math.pow(2.0,-0.5) * a;
    for(var i = 0; i < n-1; i++)
    {
        a = b;
        b = c;
        c = Math.pow(2.0,0.5) * Math.pow(i+3.0,-0.5) * x * b - Math.pow(i+2.0,0.5) * Math.pow(i+3.0,-0.5) * a;
    }
    return (n <= 0
    ?
    (-1.0) * x * Math.exp(-0.5 * x * x) * Math.pow(3.141592653589793238462,-0.25)
    :
    Math.pow(n,0.5)*Math.pow(2.0,-0.5)*a-Math.pow(n+1.0,0.5)*Math.pow(2.0,-0.5)*c);
}

function boltzQOld()
{
    let h = 0.33*0.25;
    let Ntheta = 20*4;
    let invNtheta = Math.pow(Ntheta,-1.0);
    let val = 1.0;
    let ir = this.thread.x;
    let jr = this.thread.y;
    let kr = this.thread.z;
    let hardSphere = false;
    for(let i = 0; i < 566899520; i++)
    {
        let n = i%80;
        let m = ((i-n)/80)%242;
        let l = ((((i-n)/80)-m)/242)%242;
        let k = (((((i-n)/80)-m)/242)-l)/242;
        let r = h*k;
        let vsx = h*l - 10.0;
        let vsy = h*m - 10.0;
        let theta = (2.0 * 3.141592653589793238462 * n) * invNtheta;
        let tau = Math.cos(theta)*(r-vsx) + Math.sin(theta)*(-vsy);
        let abs_tau = tau > 0.0 ? tau : -tau;
        let vpx = r - Math.cos(theta)*tau;
        let vpy = -Math.sin(theta)*tau;
        let vspx = vsx + Math.cos(theta)*tau;
        let vspy = vsy + Math.sin(theta)*tau;
        let kernel = (hardSphere ? abs_tau : 1.0);
        let prefactor = (r < 0.5*h || r > 10.0-h ? 0.5 : 1.0);
        val += prefactor * kernel * hermiteEval(Math.sqrt(vpx*vpx+vpy*vpy),2*ir) 
            * hermiteEval(Math.sqrt(vspx*vspx+vspy*vspy),2*jr) 
            * hermiteEval(r,2*kr) * 2.0 * Math.sqrt(2.0) 
            * h*h*h * invNtheta;
        val -= prefactor * kernel * hermiteEval(r,2*ir) 
            * hermiteEval(Math.sqrt(vsx*vsx+vsy*vsy),2*jr) 
            * hermiteEval(r,2*kr) * 2.0 * Math.sqrt(2.0) 
            * h*h*h * invNtheta;
    }
    return val;
}

// ir: 11, jr: 11, kr: 11, n: 80, k: 121: l: 242, m: 242
function boltzQ(ir,jr,kr,n)
{
    let h = 0.33*0.25;
    let Ntheta = 20*4;
    let invNtheta = Math.pow(Ntheta,-1.0);
    let val = 0.0;
    let k = this.thread.x;
    let l = this.thread.y;
    let m = this.thread.z;
    let hardSphere = false;
    let r = h*k;
    let vsx = h*l - 10.0;
    let vsy = h*m - 10.0;
    let theta = (2.0 * 3.141592653589793238462 * n) * invNtheta;
    let tau = Math.cos(theta)*(r-vsx) + Math.sin(theta)*(-vsy);
    let abs_tau = tau > 0.0 ? tau : -tau;
    let vpx = r - Math.cos(theta)*tau;
    let vpy = -Math.sin(theta)*tau;
    let vspx = vsx + Math.cos(theta)*tau;
    let vspy = vsy + Math.sin(theta)*tau;
    let kernel = (hardSphere ? abs_tau : 1.0);
    let prefactor = (r < 0.5*h || r > 10.0-h ? 0.5 : 1.0);
    val += prefactor * kernel * hermiteEval(Math.sqrt(vpx*vpx+vpy*vpy),2*ir) 
        * hermiteEval(Math.sqrt(vspx*vspx+vspy*vspy),2*jr) 
        * hermiteEval(r,2*kr) * 2.0 * Math.sqrt(2.0) 
        * h*h*h * invNtheta;
    val -= prefactor * kernel * hermiteEval(r,2*ir) 
        * hermiteEval(Math.sqrt(vsx*vsx+vsy*vsy),2*jr) 
        * hermiteEval(r,2*kr) * 2.0 * Math.sqrt(2.0) 
        * h*h*h * invNtheta;
    return val;
}