
package beast.evolution.util;


import test.beast.beast2vs1.trace.NumberFormatter;

import java.util.Calendar;
import java.util.TimeZone;

/**
 * A data class.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: Date.java,v 1.26 2005/05/24 20:25:57 rambaut Exp $
 */
@SuppressWarnings("Not used?")
public class Date extends TimeScale {

    public static final String DATE = "date";
    private double precision = 0.0;

    /**
     * Constructor for relative to origin
     *
     * @param time      the time in units relative to origin
     * @param units     the units of the given time
     * @param backwards true if the time is earlier than the origin
     * @param origin    the absolute origin at a Date.
     */
    public Date(double time, Type units, boolean backwards, java.util.Date origin) {
        super(units, backwards, origin);
        this.time = time;
    }

    /**
     * Constructor for an absolute date.
     *
     * @param date the date
     */
    public Date(java.util.Date date) {
        super(Units.Type.YEAR, false);
        origin = -1970.0;
        initUsingDate(date);
    }

    /**
     * Constructor for an absolute date.
     *
     * @param date  the date
     * @param units the units
     */
    public Date(java.util.Date date, Type units) {
        super(units, false);
        initUsingDate(date);
    }

    /**
     * Constructor for an absolute date with origin specified
     *
     * @param date   the date
     * @param units  the units
     * @param origin the origin as a date
     */
    public Date(java.util.Date date, Type units, java.util.Date origin) {
        super(units, false, origin);
        initUsingDate(date);
    }

    /**
     * Constructor of a relative age
     *
     * @param time      the time relative to arbitrary zero point.
     * @param units     the units the time is measured in
     * @param backwards true of the time is earlier than the zero point.
     */
    public Date(double time, Type units, boolean backwards) {
        super(units, backwards);
        this.time = time;
    }

    /**
     * Constructor for time a relative to origin
     *
     * @param origin the origin in given units from Jan 1st 1970
     */
    private Date(double time, Type units, boolean backwards, double origin) {
        super(units, backwards, origin);
        this.time = time;
    }

    //************************************************************************
    // Factory methods
    //************************************************************************

    /**
     * Create an age representing the given age (time ago) in the given units
     */
    public static Date createRelativeAge(double age, Type units) {
        return new Date(age, units, true);
    }

    /**
     * Create an age representing the given age (time ago) in the given units
     * with an origin of the given date.
     * The age represents the number units back in time from the origin.
     */
    public static Date createTimeAgoFromOrigin(double age, Type units, java.util.Date origin) {
        return new Date(age, units, true, origin);
    }

    /**
     * Create an age representing the given age (time ago) in the given units
     * with an origin as the given number of units since 1970.
     * The age represents the number units back in time from the origin.
     */
    public static Date createTimeAgoFromOrigin(double age, Type units, double origin) {
        return new Date(age, units, true, origin);
    }

    /**
     * Create an age representing the given age (time since) in the given units
     * with an origin of the given date.
     * The age represents the number units back in time from the origin.
     */
    public static Date createTimeSinceOrigin(double age, Type units, java.util.Date origin) {
        return new Date(age, units, false, origin);
    }

    /**
     * Create an age representing the given age (time since) in the given units
     * with an origin as the given number of units since 1970.
     * The age represents the number units forwards in time from the origin.
     */
    public static Date createTimeSinceOrigin(double age, Type units, double origin) {
        return new Date(age, units, false, origin);
    }

    /**
     * Create a date an sets units to Units.YEARS
     */
    public static Date createDate(java.util.Date date) {
        return new Date(date, Units.Type.YEAR);
    }

    //************************************************************************
    // Private methods
    //************************************************************************

    private void initUsingDate(java.util.Date date) {

        // get the number of milliseconds this date is after the 1st January 1970
        long millisAhead = date.getTime();


        double daysAhead = ((double) millisAhead) / MILLIS_PER_DAY;

        switch (units) {
            case DAY:
                time = daysAhead;
                break;
            case MONTH:
                time = daysAhead / DAYS_PER_MONTH;
                break;
            case YEAR:
                //time = daysAhead / DAYS_PER_YEAR;
                // more precise (so 1st Jan 2013 is 2013.0)

                // to avoid timezone specific differences in date calculations, all dates and calendars are
                // set to GMT.
                Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("GMT"));
                cal.setTime(date);

                int year = cal.get(Calendar.YEAR);
                long millis1 = cal.getTimeInMillis();

                cal.set(year, Calendar.JANUARY, 1, 0, 0);
                long millis2 = cal.getTimeInMillis();

                cal.set(year + 1, Calendar.JANUARY, 1, 0, 0);
                long millis3 = cal.getTimeInMillis();
                double fractionalYear = ((double) (millis1 - millis2)) / (millis3 - millis2);

                time = fractionalYear + year - 1970;
                break;
            default:
                throw new IllegalArgumentException();
        }


        if (time < getOrigin()) {
            time = getOrigin() - time;
            backwards = true;
        } else {
            time = time - getOrigin();
            backwards = false;
        }


    }

    /**
     * Returns the time value that is relative to the origin
     */
    public double getTimeValue() {
        return time;
    }

    /**
     * Returns the absolute time value (i.e., relative to zero).
     */
    public double getAbsoluteTimeValue() {
        if (isBackwards()) {
            return getOrigin() - getTimeValue();
        }
        return getOrigin() + getTimeValue();
    }

    public boolean before(Date date) {
        double newTime = convertTime(date.getTimeValue(), date);
        if (isBackwards()) {
            return getTimeValue() > newTime;
        }
        return getTimeValue() < newTime;
    }

    public boolean after(Date date) {
        double newTime = convertTime(date.getTimeValue(), date);
        if (isBackwards()) {
            return getTimeValue() < newTime;
        }
        return getTimeValue() > newTime;
    }

    public boolean equals(Date date) {
        double newTime = convertTime(date.getTimeValue(), date);
        return getTimeValue() == newTime;
    }

    public String getAttributeName() {
        return DATE;
    }

    public Object getAttributeValue() {
        return this;
    }

    public String toString() {
        if (isBackwards()) {
            return formatter.format(time).trim() + " " + unitString(time) + " ago";
        } else {
            return formatter.format(time).trim() + " " + unitString(time);
        }
    }

    private double time;

    private NumberFormatter formatter = new NumberFormatter(5);

    public void setPrecision(double precision) {
        this.precision = precision;
    }

    public double getPrecision() {
        return precision;
    }
}
