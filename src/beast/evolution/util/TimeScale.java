
package beast.evolution.util;

/**
 * Class for defining and converting between time scales
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TimeScale implements Units {

	/**
	 * Constructor for a timescale that starts at Jan 1st 1970. This constructor
	 * can be used if we are not interested in absolute time.
	 * param units The units
	 * param backwards True if the timescale goes backwards in time
	 */
	public TimeScale(Type units, boolean backwards) {
		this(units, backwards, 0.0);
	}
	
	/**
	 * Constructor for a timescale that starts at origin
	 * param units The units
	 * param backwards True if the timescale goes backwards in time
	 * param origin The origin specified relative to 1970 in same units
	 */
	public TimeScale(Type units, boolean backwards, double origin) {
		this.units = units;
		this.backwards = backwards;
		this.origin = origin;
	}
	
	/**
	 * Constructor for a timescale that starts at origin
	 * param units The units
	 * param backwards True if the timescale goes backwards in time
	 * param origin The origin specified as a date
	 */
	public TimeScale(Type units, boolean backwards, java.util.Date origin) {
		this.units = units;
		this.backwards = backwards;
		
		long millisAhead = origin.getTime();
		
		double daysAhead = ((double)millisAhead)/MILLIS_PER_DAY;
		
		switch (units) {	
			case DAY: this.origin = daysAhead; break;
			case MONTH: this.origin = daysAhead / DAYS_PER_MONTH; break;
			case YEAR: this.origin = daysAhead / DAYS_PER_YEAR; break;
			default: throw new IllegalArgumentException();
		}
		
	}
	
	/**
	 * @return the units of this timescale
	 */
	public Type getUnits() { return units; }
	
	/**
	 * Sets the units for this timescale.
	 */
	public void setUnits(Type units) { this.units = units; }
	
	/**
	 * @return true if larger numbers represent older dates in this time scale.
	 */
	public boolean isBackwards() { return backwards; }
	
	/**
	 * @return the date corresponding to the value zero in this time scale.
	 */
	public double getOrigin() { return origin; }
	
	/** 
	 * @return a time given in timescale scale as a time in this timescale
	 */
	public double convertTime(double time, TimeScale timeScale) {
		
		// make it forwards
		if (timeScale.isBackwards()) time = -time;
		
		// make it absolute
		time += timeScale.getOrigin();
		
		// convert to the new timescale units
		double newTime = convertTimeUnits(time, getUnits(), timeScale.getUnits());
		
		// make it relative
		newTime -= origin;
		
		// make it backwards if required
		if (backwards) newTime = -newTime;

		return newTime;
	}
	
	public String toString() {
	
		StringBuffer buffer = new StringBuffer("timescale(");
		buffer.append(unitString(0.0));
		if (backwards) {
			buffer.append(", backwards");
		} else {
			buffer.append(", forewards");
		}
		buffer.append(" from " + origin + ")");
		
		return buffer.toString();
	}
	
	public String unitString(double time) {
		String unitString = null;
		switch (units) {	
			case DAY: unitString = "day"; break;
			case MONTH: unitString = "month"; break;
			case YEAR: unitString = "year"; break;
			default: throw new IllegalArgumentException();
		}
		if (time == 1.0) {
			return unitString;
		} else return unitString + "s";
	}
	
	public static void main(String[] args) {
	
		TimeScale timeScale1 = new TimeScale(Type.DAY, true);
		TimeScale timeScale2 = new TimeScale(Type.YEAR, true);
		
		System.out.println(timeScale1);
		System.out.println(timeScale2);
		
		double testTime = 100.0;
		System.out.println("Test time = " + testTime);
		
		System.out.println("timeScale1.convertTime(" + testTime + ", timeScale2)=" + timeScale1.convertTime(testTime, timeScale2));
		System.out.println("timeScale2.convertTime(" + testTime + ", timeScale1)=" + timeScale2.convertTime(testTime, timeScale1));
	}
	
	
	//*************************************************************************
	// STATIC STUFF
	//*************************************************************************
	
	/** 
	 * @return time in currentUnits as newUnits.
	 */
	public static double convertTimeUnits(double time, Type currentUnits, Type newUnits) {
		
		return time * getScale(currentUnits, newUnits);
	}
	
	/** 
	 * @return the scaling factor for converting currentUnits into newUnits.
	 */
	public static double getScale(Type currentUnits, Type newUnits) {
		if (currentUnits == newUnits) return 1.0;
		
		switch (currentUnits) {
			case DAY:
				switch (newUnits) {
					case MONTH: return 1.0/DAYS_PER_MONTH;
					case YEAR: return 1.0/DAYS_PER_YEAR;
					default: throw new IllegalArgumentException();
				}
			case MONTH:
				switch (newUnits) {
					case DAY: return DAYS_PER_MONTH;
					case YEAR: return 1.0/MONTHS_PER_YEAR;
					default: throw new IllegalArgumentException();
				}
			case YEAR:
				switch (newUnits) {
					case DAY: return DAYS_PER_YEAR;
					case MONTH: return MONTHS_PER_YEAR;
					default: throw new IllegalArgumentException();
				}
			default: throw new IllegalArgumentException();
		}
	}
	
	//*************************************************************************
	// PRIVATE STUFF
	//*************************************************************************

	// The origin is specified in days relative to 1st January 1970
	protected double origin = 720035.0;
	protected Type units;
	protected boolean backwards;
	
	protected static double MILLIS_PER_DAY = 86400000.0;
	protected static double DAYS_PER_YEAR = 365.25;
	protected static double MONTHS_PER_YEAR = 12.0;
	protected static double DAYS_PER_MONTH = DAYS_PER_YEAR / MONTHS_PER_YEAR;
}
