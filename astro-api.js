//file courtesy of http://www.astro-phys.com/js/astro/api.js
//it had a private function (apiCall) that I had to modify

var astro = (function(astro) {
  astro.JSURL = 'http://www.astro-phys.com/js/astro/'
  astro.URL = 'http://www.astro-phys.com/api/';

  function apiCall(type, args, callback) {
	 $.ajax({
		dataType: "jsonp",
		url: astro.URL + type + '?callback=?',
		data: args,
		success: callback 
	});
  }
  astro.apiCall = apiCall;

  astro.constants = {denum:406,lenum:406,tdatef:0,tdateb:11997060422593800,center:0,clight:299792.458,au:1.49597870691E8,emrat:81.30056,gm1:4.912547451450812E-11,gm2:7.243452486162703E-10,gmb:8.997011346712499E-10,gm4:9.549535105779258E-11,gm5:2.8253459095242264E-7,gm6:8.459715185680659E-8,gm7:1.2920249167819694E-8,gm8:1.5243589007842763E-8,gm9:2.1886997654259697E-12,gms:2.959122082855911E-4,rad1:2439.76,rad2:6052.3,rad4:3397.515,jdepoc:2440400.5,x1:0.3617627146035093,y1:-0.09078196772958605,z1:-0.08571498318176335,
xd1:0.00336749391398414,yd1:0.024894520467648874,zd1:0.012946300688650358,x2:0.612751941341838,y2:-0.3483653684949682,z2:-0.19527828898022992,xd2:0.010952068361699194,yd2:0.015617684365262055,zd2:0.006331105553601057,xb:0.12051741729540039,yb:-0.9258384789394528,zb:-0.401540220801811,xdb:0.016811268303992644,ydb:0.0017483093133816164,zdb:7.582028689905378E-4,x4:-0.11018607428285881,y4:-1.3275994561325557,z4:-0.6058891326142037,xd4:0.014481653059735103,yd4:2.424631177602962E-4,zd4:-2.8152073424800533E-4,
x5:-5.379706898835912,y5:-0.8304805814601515,z5:-0.22482870022837284,xd5:0.0010920115430088524,yd5:-0.006518116565792683,zd5:-0.0028207831653650475,x6:7.894392441979052,y6:4.596477801626945,z6:1.5586975735302668,xd6:-0.003217555239303309,yd6:0.004335809858955276,zd6:0.0019286465656538445,x7:-18.265398306822032,y7:-1.1619445055181146,z7:-0.25010348393737786,xd7:2.2118841741776543E-4,yd7:-0.003762475932846577,zd7:-0.001651014703068002,x8:-16.05504258376823,y8:-23.942181216178568,z8:-9.40015672354883,
xd8:0.0026427710433669163,yd8:-0.0014983144553592097,zd8:-6.790419030179561E-4,x9:-30.483319603999316,y9:-0.8724783554957969,z9:8.911563040989932,xd9:3.2221044772358804E-4,yd9:-0.0031435703021520205,zd9:-0.0010779488297393042,xm:-8.081773279114842E-4,ym:-0.0019946300016203994,zm:-0.0010872626608381018,xdm:6.010848166591299E-4,ydm:-1.6744546061515147E-4,zdm:-8.556214497398616E-5,xs:0.0045025081562338936,ys:7.670747009323788E-4,zs:2.6605680517702713E-4,xds:-3.5174820964518867E-7,yds:5.17762539958483E-6,
zds:2.229101854391665E-6,beta:1,gamma:1,j2sun:2.0E-7,gdot:0,ma0001:1.390787378942278E-13,ma0002:2.959122082855911E-14,ma0004:3.846858707712684E-14,mad1:1.8,mad2:2.4,mad3:5,re:6378.137,asun:696E3,phi:0.0051299597051581245,tht:0.3823906558768601,psi:1.2941422241102787,omegax:4.5247044990228E-5,omegay:-2.23092763198743E-6,omegaz:0.229944858701367,am:1738,j2m:2.0431200665465293E-4,j3m:8.78546950779311E-6,j4m:-1.45383007072E-7,c22m:2.251782439166225E-5,c31m:3.0803809783429375E-5,c32m:4.87980736275513E-6,
c33m:1.7701764624348135E-6,s31m:4.259328630030294E-6,s32m:1.6955162680003687E-6,s33m:-2.7097009665669175E-7,c41m:-7.17780149806E-6,c42m:-1.43951838385E-6,c43m:-8.54788154819E-8,c44m:-1.5490389313E-7,s41m:2.94743374914E-6,s42m:-2.8843721272E-6,s43m:-7.88967312839E-7,s44m:5.6404155572E-8,lbet:6.316121342194731E-4,lgam:2.2785831399775596E-4,k2m:0.029922116659705348,taum:0.16671655584928446,ae:6378.137,j2e:0.001082626,j3e:-2.533E-6,j4e:-1.616E-6,k2e0:0.34,k2e1:0.3,k2e2:0.3,taue0:0,taue1:0.01290895939156023,
taue2:0.006941785584052303,drotex:2.44E-4,drotey:-0.001193,gmast1:6.466825433842555E-14,gmast2:1.277481189104146E-14,gmast3:3.334058772960295E-15,kvc:0,ifac:3.0E-4,phic:-0.0042595183,thtc:0.4088443,psic:-1.714509,omgcx:0,omgcy:-1.5816707E-6,omgcz:0.229888,psidot:0,mgmis:1,rotex:0,rotey:0};

  astro.julianToCalendar = function(julian) {
    var jd = 0.5 + julian;
    var I = Math.floor(jd);
    var F = jd - I;
    if (I > 2229160) {
      var A = Math.floor((I - 1867216.25) / 36524.25);
      var B = I + 1 + A - Math.floor(A / 4.0);
    } else {
      var B = I;
    }
    var C = B + 1524;
    var D = Math.floor((C - 122.1) / 365.25);
    var E = Math.floor(365.25 * D);
    var G = Math.floor((C - E) / 30.6001);
    var d = C - E + F - Math.floor(30.6001 * G);
    if (G < 13.5) {
      var month = G - 1;
    } else {
      var month = G - 13;
    }
    if (month > 2.5) {
      var year = D - 4716;
    } else {
      var year = D - 4715;
    }
    var day = Math.floor(d);
    var h = (d - day) * 24;
    var hour = Math.floor(h);
    var mind = Math.abs(h - hour) * 60;
    var min = ''+Math.floor(mind);
    if (min.length == 1) min = '0' + min;
	var sec = ''+Math.floor((mind - min) * 60);
    if (sec.length == 1) sec = '0' + sec;
	return year + '-' + month + '-' + day + ' ' + hour + ':' + min + ':' + sec;
  };

  function State(position, velocity) {
    this.position = position;
    this.velocity = velocity;
  }
  astro.State = State;

  State.prototype.scale = function(factor) {
    for(var i = 0; i < 3; ++i) {
      this.position[i] *= factor;
      this.velocity[i] *= factor;
    }
  }

  astro.getStates = function(bodies, date, callback) {
    var query = {bodies: bodies.join(','), date: date};
    apiCall('states', query, function(data) {
      for(var body in data.results) {
        var state = data.results[body];
        data.results[body] = new State(state[0], state[1]);
      }
      callback(data.results, data.date);
    });
  }

  astro.getState = function(body, date, callback) {
    astro.getStates([body], date, function(results, date) {
      callback(results[body], date);
    });
  }

  function chebyshevEvaluate(x, coeffs) {
    var x2 = 2.0 * x, d = 0.0, dd = 0.0, tmp;
    for(var i = coeffs.length - 1; i >= 1; --i) {
      tmp = d;
      d = x2 * d - dd + coeffs[i];
      dd = tmp;
    }
    return x * d - dd + coeffs[0];
  }

  function chebyshevDerive(x, coeffs) {
    var x2 = 2.0 * x, d = 0.0, dd = 0.0, tmp;
    for(var i = coeffs.length - 1; i >= 2; --i) {
      tmp = d;
      d = x2 * d - dd + k * coeffs[i];
      dd = tmp;
    }
    return x2 * d - dd + coeffs[1];
  }

  function CoeffSet(coeffs, start, end) {
    this.coeffs = coeffs;
    this.start = start;
    this.end = end;
    this.step = end - start;
  }
  astro.CoeffSet = CoeffSet;

  CoeffSet.prototype.getPosition = function(date) {
    if(date < this.start || date > this.end)
      throw new Error('date out of range for Coeffs');
    var x = 2.0 * (date - this.start) / this.step - 1.0;
    var third = this.coeffs.length / 3;
    var position = [];
    for(var i = 0; i < 3; ++i) {
      var values = this.coeffs.slice(i * third, (i + 1) * third);
      position[i] = chebyshevEvaluate(x, values);
    }
    return position;
  }

  CoeffSet.prototype.getState = function(date) {
    if(date < this.start || date > this.end)
      throw new Error('date out of range for Coeffs');
    var x = 2.0 * (date - this.start) / this.step - 1.0;
    var step2 = 2.0 / this.step;
    var third = this.coeffs.length / 3;
    var position = [], velocity = [];
    for(var i = 0; i < 3; ++i) {
      var values = this.coeffs.slice(i * third, (i + 1) * third);
      position[i] = chebyshevEvaluate(x, values);
      velocity[i] = chebyshevDerive(x, values) * step2;
    }
    return new astro.State(position, velocity);
  }

  astro.getCoeffSets = function(bodies, date, callback) {
    var query = {bodies: bodies.join(','), date: date};
    apiCall('coeffs', query, function(data) {
      for(var body in data.results) {
        var coeffs = data.results[body].coeffs;
        var start = data.results[body].start;
        var end = data.results[body].end;
        data.results[body] = new CoeffSet(coeffs, start, end);
      }
      callback(data.results, data.date);
    });
  }
  
  astro.getCoeffSet = function(body, date, callback) {
    astro.getCoeffSets([body], date, function(results, date) {
      callback(results[body], date);
    });
  }

  function Record(coeffs, start, end) {
    this.coeffs = coeffs;
    this.start = start;
    this.end = end;
    this.step = this.end - this.start;
  }
  astro.Record = Record;

  Record.prototype.getPosition = function(body, date) {
    if(date < this.start || date > this.end)
      throw new Error('date out of range for Coeffs');
    var coeffs = this.coeffs[body];
    var step = (this.end - this.start);
    var index = parseInt(((date - this.start) / step) * coeffs.length);
    return coeffs[index].getPosition(date);
  }

  Record.prototype.getPositions = function(date) {
    if(date < this.start || date > this.end)
      throw new Error('date out of range for Coeffs');
    var results = {};
    var earthmoon = this.getPosition('earthmoon', date);
    var geomoon = this.getPosition('geomoon', date);
    var emrat = 1.0 / (1.0 + astro.constants.emrat);
    var earth = [earthmoon[0] - geomoon[0] * emrat,
                 earthmoon[1] - geomoon[1] * emrat,
                 earthmoon[2] - geomoon[2] * emrat];
    var moon = [geomoon[0] + earth[0], geomoon[1] + earth[1], geomoon[2] + earth[2]];
    results['sun']     = this.getPosition('sun', date);
    results['mercury'] = this.getPosition('mercury', date);
    results['venus']   = this.getPosition('venus', date);
    results['moon']    = moon;
    results['earth']   = earth;
    results['mars']    = this.getPosition('mars', date);
    results['jupiter'] = this.getPosition('jupiter', date);
    results['saturn']  = this.getPosition('saturn', date);
    results['uranus']  = this.getPosition('uranus', date);
    results['neptune'] = this.getPosition('neptune', date);
    results['pluto']   = this.getPosition('pluto', date);
    return results;
  }

  Record.prototype.getStates = function(date, callback) {
    if(date < this.start || date > this.end)
      throw new Error('date out of range for Coeffs');
    for(var body in this.coeffs) {
      var coeffs = this.coeffs[body];
      var step = (this.end - this.start) / coeffs.length;
      var index = parseInt((date - this.start) / step);
      callback(body, coeffs[index].getState(date));
    }
  }

  astro.getRecord = function(date, callback) {
    apiCall('records', {date:date}, function(data) {
      var results = data.results;
      var start = data.start;
      var end = data.end;
      for(var body in results) {
        var nchunks = results[body].length
        var c_step = (end - start) / nchunks;
        for(var i = 0; i < nchunks; ++i) {
          var c_start = start + i * c_step;
          var c_end = c_start + c_step;
          var coeffs = results[body][i];
          results[body][i] = new astro.CoeffSet(coeffs, c_start, c_end);
        }
      }
      callback(new Record(results, start, end), data.date);
    });
  }
  
  astro.now = function() {
    var date = new Date();
    var year = date.getFullYear();
    var month = date.getMonth() + 1;
    var day = date.getDate();
    var hour = date.getHours();
    var min = date.getMinutes();
    var sec = date.getSeconds();
    return year + '-' + month + '-' + day + ' ' + hour + ':' + min + ':' + sec;
  };
  

  return astro;
}(astro || {}));


